from typing import List, Optional
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from tracetrack.entities.record import TraceSeqRecord, Record, Feature
from tracetrack.entities.aligned import AlignedReference, AlignedTrace
from tracetrack.entities.sequence_position import SequencePosition, UNKNOWN_INSERTION
from tracetrack.alignment_utils import align_clustalo, suffix

MUTATION_TYPE_NONE = 'no_mut'
MUTATION_TYPE_SILENT = 'silent_mut'
MUTATION_TYPE_MISSENSE = 'mis_mut'
MUTATION_TYPE_NONSENSE = 'nons_mut'
MUTATION_TYPE_FRAMESHIFT = "frameshift"


class Alignment:
    """
    Class for storing results of an alignment. Contains both reference and aligned sequences, mutations and additional
    metrics like coverage.
    """
    def __init__(self, aligned_reference: AlignedReference, aligned_traces: List[AlignedTrace], population: str = None):
        """
        :param aligned_reference: aligned reference record
        :param aligned_traces: list of aligned trace records
        :param population: Name of aligned population
        """
        self.aligned_reference = aligned_reference
        self.aligned_traces = aligned_traces
        self.population = population
        self.cds = [feat.location for feat in aligned_reference.record.features if feat.type == "CDS"]
        self.positions = self._create_sequence_positions()

        #
        # Pre-calculate statistics:
        #

        self.covered_positions = sum(1 for pos in self.positions if pos.has_coverage())
        self.mutations = sum(1 for pos in self.positions if pos.is_mutated())
        # SequencePosition can't know about frameshift

        mutation_list = [self.mutation_type(pos) for pos in range(len(self.positions))]
        self.mismatches = [pos for pos in range(len(self.positions)) if self.mutation_type(pos) != MUTATION_TYPE_NONE]
        self.silent_mutations = sum(t == MUTATION_TYPE_SILENT for t in mutation_list)
        self.mis_mutations = sum(t == MUTATION_TYPE_MISSENSE for t in mutation_list)
        self.nons_mutations = sum(t == MUTATION_TYPE_NONSENSE for t in mutation_list)
        self.frameshifts = sum(t == MUTATION_TYPE_FRAMESHIFT for t in mutation_list)
        self.mutations = self.silent_mutations + self.mis_mutations + self.frameshifts + self.nons_mutations

        self.perc_coverage = (self.covered_positions / len(self.positions)) * 100
        self.perc_identity = \
            ((self.covered_positions - self.silent_mutations - self.mis_mutations - self.nons_mutations)
             / self.covered_positions) * 100

        self.kozak_count = sum(1 for aligned_trace in self.aligned_traces if aligned_trace.has_kozak())
        self.stop_count = sum(1 for aligned_trace in self.aligned_traces if aligned_trace.has_stop())

    @classmethod
    def align_multiple(cls, ref_record: Record, trace_records: List[TraceSeqRecord], population: str = None) -> 'Alignment':
        """
        Align all sequences in list to a reference seq
        :param ref_record: reference record to align
        :param trace_records: list of trace sequence records to align
        :param population: Name of aligned population
        :param cds: Optional list of CDS features from GenBank file
        :return: Alignment
        """
        ref_aligned_seq, query_aligned_seqs = align_clustalo(ref_record, trace_records)

        # get reference start and end without gaps
        ref_start = 0
        while ref_aligned_seq[ref_start] == '-':
            ref_start += 1
        ref_end = len(ref_aligned_seq)
        while ref_end > 0 and ref_aligned_seq[ref_end - 1] == '-':
            ref_end -= 1

        aligned_reference = AlignedReference(ref_record, ref_aligned_seq[ref_start:ref_end])

        # get kozak and stop sequence from ends
        aligned_traces = []
        for trace_record, query_align in zip(trace_records, query_aligned_seqs):
            leading_seq = str(query_align[:ref_start]).replace('-', '')
            trailing_seq = str(query_align[ref_end:]).replace('-', '')

            # TODO add test for leading and trailing region with multiple traces
            aligned_traces.append(AlignedTrace(
                record=trace_record,
                aligned_seq=query_align[ref_start:ref_end],
                leading_seq=leading_seq,
                trailing_seq=trailing_seq
            ))

        return cls(
            aligned_reference=aligned_reference,
            aligned_traces=aligned_traces,
            population=population
        )

    @property
    def ref_id(self):
        return self.aligned_reference.record.id

    def _create_sequence_positions(self) -> List['SequencePosition']:
        """
        For each position in reference sequence, create list of possible options from reads and create
        SequencePosition object. Skip dashes and N's.
        A conflict happens when there is a double peak in all trace sequences / in the only trace sequence
        :return: list of sequence position objects
        """
        seqlist = []
        result_bases = 0
        ref_bases = 0
        for i, ref_base in enumerate(self.aligned_reference.aligned_seq):
            relevant_traces = [t for t in self.aligned_traces if t.is_considered(i)]
            num_mixed_peaks = sum(i in t.mixed_peaks for t in relevant_traces)
            mixed = num_mixed_peaks >= 2 or num_mixed_peaks == len(relevant_traces) == 1
            seqlist.append(SequencePosition(
                options=[t.aligned_seq[i] for t in relevant_traces],
                ref=ref_base,
                pos_ref=ref_bases,
                pos_mutated=result_bases,
                mixed_peak=mixed
            ))
            if ref_base != "-":
                ref_bases += 1
            if seqlist[-1].result != "-" and seqlist[-1].result != "?":
                # unknown insert is not counted as a valid position
                result_bases += 1
        return seqlist

    def is_gap(self, pos: int, mutated: bool) -> bool:
        """
        Checks if there is a gap  at given position in mutated or reference sequence.
        :param pos: int
        :param mutated: bool
        :return: bool
        """
        if mutated:
            return self.positions[pos].result == "-"
        else:
            return self.positions[pos].ref == "-"

    def get_zero_coverage_runs(self):
        """
        Finds all regions of alignment with zero coverage
        :return: List of locations (tuple of start and end) in one-based counting
        """
        zeros = []
        start = 0
        running = False
        for i, pos in enumerate(self.positions):
            if pos.has_coverage():
                if running:
                    zeros.append((start+1, i))
                    # numbering starts from 1
                running = False
            else:
                if not running:
                    start = i
                running = True
        if running:
            zeros.append((start+1, len(self.positions)))
        return zeros

    def get_codon_reference(self, pos: int) -> Optional[Seq]:
        """Get the codon at specfied position in the reference sequence."""
        return self.get_codon(pos, False)

    def get_codon_mutated(self, pos: int) -> Optional[Seq]:
        """Get the codon at specfied position in the mutated sequence."""
        return self.get_codon(pos, True)

    def get_codon(self, pos: int, mutated: bool) -> Optional[Seq]:
        """
        Return string of three bases belonging to codon at given sequence position.
        :param pos: int
        :param mutated: bool
        :return: Seq
        """
        num = self.get_codon_num(pos, mutated)
        codon = []
        left = pos
        while left >= 0 and self.get_codon_num(left, mutated) == num:
            if not self.is_gap(left, mutated):
                codon.insert(0, self.get_base_at_position(left, mutated))
            left -= 1
        right = pos + 1
        while right < len(self.positions) and self.get_codon_num(right, mutated) == num:
            if not self.is_gap(right, mutated):
                codon.append(self.get_base_at_position(right, mutated))
            right += 1
        return Seq("".join(codon))

    # might be a place to use lru_cache?
    def get_base_at_position(self, pos: int, mutated: bool) -> str:
        if mutated:
            return self.positions[pos].result
        else:
            return self.positions[pos].ref

    # might be a place to use lru_cache?
    def is_same_aminoacid(self, pos: int) -> bool:
        """
        Return bool value specifying whether the reference and mutated sequence code for the same amino acid
        at the given position.
        """
        try:
            AA1 = self.translate(pos, False)
            AA2 = self.translate(pos, True)
        except TranslationError:
            return False
        # If there are not enough bases to translate, AA2 is "" and therefore False is returned
        return AA1 == AA2

    def translate_reference(self, pos: int) -> str:
        """Translate the DNA codon at the given position in the reference sequence to an amino acid."""
        return self.translate(pos, False)

    def translate_mutated(self, pos: int) -> str:
        """Translate the DNA codon at the given position in the mutated sequence to an amino acid."""
        return self.translate(pos, True)

    def translate(self, pos: int, mutated: bool) -> str:
        """
        Translate the DNA codon at the given position in the sequence to an amino acid. Ignore question marks in the
        sequence (unknown insertions).
        """
        codon = self.get_codon(pos, mutated)
        if codon is None:
            return ""
        if UNKNOWN_INSERTION in str(codon):
            codon = Seq(str(codon).replace("?", ""))
        # if the codon is empty (reached end of sequence) or just < 3 bases remaining, AA is ""
        AA = codon.translate()
        return AA

    # might be a place to use lru_cache?
    def mutation_type(self, pos: int) -> str:
        """
        Compare bases at given position, return correct option for type of mutation. If the bases differ, is_same_aminoacid
        is called to determine mutation type.
        :param pos: int
        :return: string (one of options defined at beginning of file)
        """
        position = self.positions[pos]
        ref = position.ref.upper()
        ch = position.result.upper()

        feat = self.aligned_reference.get_feature_for_pos(position.pos_ref)

        if not feat.is_coding():
            return MUTATION_TYPE_NONE if ch == ref else MUTATION_TYPE_SILENT

        else:
            ref_offset, mut_offset = self.get_offset(feat.location.start)
            ref_pos = position.pos_ref - ref_offset
            mut_pos = position.pos_mutated - mut_offset
            if (ref_pos - mut_pos) % 3 != 0 and position.result != "-" and position.coverage > 0:
                return MUTATION_TYPE_FRAMESHIFT
            elif pos + 1 < len(self.positions):
                next_pos = self.positions[pos + 1]
                if position.is_indel() and (next_pos.pos_ref - ref_offset - next_pos.pos_mutated + mut_offset) % 3 != 0:
                    return MUTATION_TYPE_FRAMESHIFT
        if ch == ref:
            return MUTATION_TYPE_NONE
        if self.is_same_aminoacid(pos):
            return MUTATION_TYPE_SILENT
        try:
            translation = self.translate(pos, True)
        except TranslationError:
            return MUTATION_TYPE_NONSENSE
        if translation == '*':
            return MUTATION_TYPE_NONSENSE
        return MUTATION_TYPE_MISSENSE

    def get_trace_ids(self):
        return [aligned_trace.record.id for aligned_trace in self.aligned_traces]

    def format_position(self, pos):
        """
        Return a format string which is used in the html output for a given base.
        """
        feat = self.aligned_reference.get_feature_for_pos(self.positions[pos].pos_ref)
        format_str = feat.get_format()

        mutation = self.mutation_type(pos)
        if self.positions[pos].mixed_peak:
            mutation = mutation + " mixed_peak"
        return format_str + " " + mutation

    def get_codon_num(self, pos: int, mutated: bool) -> int:
        """
        Get the position number of the codon at the given sequence position.
        """
        ref_pos = self.positions[pos].pos_ref
        feat = self.aligned_reference.get_feature_for_pos(ref_pos)
        ref_offset, mut_offset = self.get_offset(feat.location.start)
        if mutated:
            return (self.positions[pos].pos_mutated - mut_offset) // 3
        return (self.positions[pos].pos_ref - ref_offset) // 3

    def get_offset(self, start):
        """Get the position of the start of coding sequence in the alignment."""
        cds_start_in_mut = self.aligned_reference.ref_to_mut[start]
        mut_offset = self.positions[cds_start_in_mut].pos_mutated
        return start, mut_offset

    def position_to_display(self, pos):
        """
        Produce a position label in the reference sequence from a position in the mutated sequence to be displayed in tooltip.
        For insertions, insertion code is used by adding letters to the previous reference position. E.g. if two bases
        are inserted after position 11, they will be labeled 11a and 11b.
        For all other positions, the label is the same as the reference position it aligns to.
        :param pos: position in alignment
        :return: position in the reference sequence as string
        """
        i = pos
        # find last position aligned to reference without insertion
        while self.positions[i].ref in ["-", UNKNOWN_INSERTION]:
            i -= 1
        # number of inserted positions
        diff = pos - i
        letter_code = suffix(diff)
        # return number along with letter code
        return str(self.positions[i].pos_ref + 1) + letter_code

    def get_translation(self):
        """
        Produce complete translation to amino acids from alignment. Non-coding regions are translated to spaces.
        """
        codon_numbers_ref = [self.get_codon_num(i, False) for i in range(len(self.positions))]
        coding = [self.aligned_reference.get_feature_for_pos(i).is_coding() for i in range(len(self.positions))]
        translation = []
        for pos in range(len(self.positions)):
            if not coding[pos]:
                translation.append(" ")
            elif self.positions[pos].ref == "-":
                continue
            else:
                if (pos == 0) or (codon_numbers_ref[pos - 1] < codon_numbers_ref[pos]) or (translation[-1] == "┘") or (not coding[pos - 1]):
                    # first nt of codon
                    translation.append("└")
                elif (pos == len(self.positions) - 1) or (codon_numbers_ref[pos + 1] > codon_numbers_ref[pos]):
                    # last nt of codon
                    translation.append("┘")
                else:
                    ref_aa = self.translate_reference(pos)
                    mut_aa = self.translate_mutated(pos)
                    if ref_aa == mut_aa:
                        translation.append(str(ref_aa))
                    else:
                        translation.append(f"{ref_aa} → {mut_aa}")
        return translation

