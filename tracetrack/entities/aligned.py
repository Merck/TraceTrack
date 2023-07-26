from tracetrack.entities.record import Record, TraceSeqRecord, Feature
from tracetrack.alignment_utils import preserve_lowercase, translate_codon
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from typing import List


NOT_CONSIDERED_BASES = ['-', 'N']
KOZAK = "GCCACC"


class AlignedReference:
    """Class for storing the reference sequence used in an alignment."""
    def __init__(self, record: Record, aligned_seq: Seq):
        self.record = record
        self.aligned_seq, self.ref_to_mut = preserve_lowercase(record.seq, aligned_seq)

    def get_feature_for_pos(self, pos):
        """
        Return the feature present at given position (e.g. coding sequence) or a default feature when none is present.
        """
        possible = [feat for feat in self.record.features if feat.location.start <= pos < feat.location.end]
        if len(possible) == 0:
            # TODO: should each position be part of a feature?
            return Feature("default", FeatureLocation(pos, pos+1, 1), "none", -1)
        cds = [feat for feat in possible if feat.type == "CDS"]
        if len(cds) == 0:
            return possible[0]
        if len(cds) > 1:
            pass
            # TODO: raise some error when CDS locations overlap?
        return cds[0]


class AlignedTrace:
    def __init__(self, record: TraceSeqRecord, aligned_seq: Seq, leading_seq: str, trailing_seq: str):
        """
        This alignment:
         ref: ----------CGTCCCTACACAAGTTC-----
        seq1: NNTAC-TGAGCGTCCC---ACAAGTTCT-GNN
        seq2: NNTACATGAGCGTCCC---ACAAGTTCTGGNN

        Then seq1 will be stored as:
          record.seq: "NNTACTGAGCGTCCCACAAGTTCTGNN"
         leading_seq: "NNTACTGAG"
         aligned_seq: "CGTCCC---ACAAGTTC"
        trailing_seq: "TGNN"

        :param record: trace sequence record of origin
        :param aligned_seq: aligned sequence (with gaps), without regions outside aligned reference
        :param leading_seq: 5' leading sequence (without gaps) from 0 up to 6 nucleotides (used to recognize kozak sequence)
        :param trailing_seq: 3' trailing sequence (without gaps) from 0 up to 3 nucleotides (used to recognize stop codon)
        """
        self.record = record
        self.aligned_seq = aligned_seq
        self.leading_seq = leading_seq
        self.trailing_seq = trailing_seq
        self.first_base_idx = self._get_first_base_idx()
        self.last_base_idx = self._get_last_base_idx()
        self.mixed_peaks = self.find_double_peaks()

    def has_kozak(self):
        return self.leading_seq[-6:] == KOZAK

    def find_double_peaks(self):
        """
        Get a list of positions with mixed peaks.
        """
        mapping = self.get_aligned_positions()
        doubles = []

        for pos in range(len(self.aligned_seq)):
            if mapping[pos] in self.record.mixed_peaks:
                doubles.append(pos)
        return doubles

    def has_stop(self):
        """Indicates whether the aligned sequence ends with a stop codon."""
        return len(self.trailing_seq) >= 3 and translate_codon(self.trailing_seq[:3]) == "*"

    def _get_first_base_idx(self):
        """Returns first considered base in aligned_sequence"""
        for idx, base in enumerate(self.aligned_seq):
            if base not in NOT_CONSIDERED_BASES:
                return idx
        return None

    def _get_last_base_idx(self):
        """Returns last considered base in aligned_sequence"""
        for idx_from_end, base in enumerate(self.aligned_seq[::-1]):
            if base not in NOT_CONSIDERED_BASES:
                return len(self.aligned_seq) - idx_from_end - 1
        return None

    def is_considered(self, i):
        """
        Check if given position should be considered for voting on result.
        Returns true for all valid bases and for gaps inside the sequence.
        Returns false for Ns and for gaps at edges.
        :param i:
        :return: whether given position should be considered for voting on result
        """
        # Ignore Ns
        if self.aligned_seq[i] == 'N' or self.aligned_seq[i] == 'n':
            return False

        # Only gaps in aligned part
        if self.first_base_idx is None:
            return False

        # Ignore gaps at edges
        if self.aligned_seq[i] == '-' and (i < self.first_base_idx or i > self.last_base_idx):
            return False

        return True

    def get_aligned_positions(self) -> List[int]:
        """
        Get mapping from aligned sequence position to original sequence position
        :return: List where list[alignedPos] = origPos
        """
        positions = []
        pos = len(self.leading_seq)
        for base in self.aligned_seq:
            positions.append(pos)
            if base != '-':
                pos += 1

        return positions

