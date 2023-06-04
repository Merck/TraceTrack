from tracetrack.entities.record import TraceSeqRecord, Record
from tracetrack.entities.alignment import Alignment
from tracetrack.alignment_utils import align_global_score
from tracetrack.server_utils import has_extension, ALLOWED_REF_EXTENSIONS
from tracetrack.entities.errors import AlignmentError
from openpyxl import load_workbook
from typing import List, Tuple
from flask import flash
from csv import reader
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def read_csv(filename):
    """
    Read a reference sheet in csv format
    """
    with open(filename, 'r') as f:
        csv_reader = reader(f)
        header = next(csv_reader)
        records = [Record.read_from_csv(row) for row in csv_reader]
    return records, header


def read_xlsx(filename):
    """
    Read a reference sheet in xlsx format
    """
    table = pd.read_excel(filename, dtype=str, engine='openpyxl')
    header = table.columns[:2]
    records = [Record.read_from_csv(row) for i, row in table.iterrows() if row[0] and not pd.isna(row[0])]

    return records, header


def read_fasta(filename):
    """
    Read reference sheet in fasta format
    """
    iterator = SeqIO.parse(filename, "fasta")
    records = [Record.read_from_fasta(seq_record) for seq_record in iterator]
    return records


def read_gb(filename):
    """
    Read reference sheet in GenBank formati
    """
    gb_records = SeqIO.parse(open(filename, "r"), "genbank")
    records = [Record.read_from_genbank(gb_record) for gb_record in gb_records]
    return records


class ReferenceDb:
    """
    Class for storing reference sequences loaded from reference sheet provided by user.
    Each reference sequence is stored as a Record object.
    """
    def __init__(self, records):
        self.sequences = records

        if len(records) != len(set(records)):
            flash(f"There are multiple references with same ID. The first one will be used.")
            # FIXME

    @classmethod
    def read_file(cls, filename) -> 'ReferenceDb':
        """
        Read a DataFrame from an .xlsx, .csv, fasta or .gbk file. Return ReferenceDb object.
        :param filename: str
        :return: ReferenceDb
        """
        if has_extension(filename, ["csv", "xlsx", "xls"]):
            if has_extension(filename, ["csv"]):
                records, header = read_csv(filename)
            else:
                records, header = read_xlsx(filename)
            if len(records) > 0 and header[0].lower() != 'id':
                flash('Warning: Expected first column name “ID”, got “{}”, make sure to include a valid header.'.format(
                    header[0]))
            if len(records) > 0 and header[1].lower() != 'sequence':
                flash(
                    'Warning: Expected  second column name “Sequence”, got “{}”, make sure to include a valid header.'.format(
                        header[1]))

        elif has_extension(filename, ["fa", "fasta"]):
            records = read_fasta(filename)

        elif has_extension(filename, ["genbank", "gbk", "gb"]):
            records = read_gb(filename)

        else:
            raise ReferenceError(f"The provided file {filename} has none of the accepted extensions: {ALLOWED_REF_EXTENSIONS}")
            # TODO

        return cls(records)

    def get_ref_sequence(self, id) -> Record:
        """
        Return the reference sequence record from table based on its ID. Returns None if the ID is not found.
        :param id: reference ID
        :return: reference sequence Record
        """
        rec = next((x for x in self.sequences if x.id == id), None)
        return rec

    # useless
    def create_ref_record(self, id) -> Record:
        """
        Create reference sequence record
        :param id: reference ID
        :return: SeqRecord with sequence and reference ID
        """
        return Record(Seq(self.get_ref_sequence(id).seq), id)

    def align_sequences(self, records: List[TraceSeqRecord], refid, population: str = None) -> Tuple[
        Alignment, List[str]]:
        """
        Sort list of TraceSeqRecord by reference id, create Alignment object for each reference, return as a list.
        :param records: list of TraceSeqRecords to align
        :param population: str, name of population
        :param search: bool, align each trace sequence to the best matching reference instead of parsing filename
        :return: list of Alignment objects, one object for each reference
        """
        ref_record = self.get_ref_sequence(refid)
        alignment = None
        warnings = []
        try:
            alignment = Alignment.align_multiple(ref_record, records, population=population)
        except AlignmentError as e:
            warnings.append(f"Alignment error for reference ID {refid}: {e.message}")

        return alignment, warnings

    def match_trace_to_ref(self, record: TraceSeqRecord):
        """
        Assign given TraceSeqRecord a reference id and direction by best match.
        :param record: Trace record to find the best reference for
        :return: nothing
        """
        # match according to filename
        max_match = 0
        reference = None
        for ref in self.sequences:
            if ref.id in record.id:
                if len(ref.id) > max_match:
                    max_match = len(ref.id)
                    reference = ref
        if reference is not None:
            fwd = align_global_score(reference.seq, record.seq)
            rev = align_global_score(reference.seq, record.seq.reverse_complement())
            if rev > fwd:
                record.flag_as_reverse(True)
            record.assign_reference(reference.id)

        # match by best alignment score
        else:
            scores_fwd = [align_global_score(ref.seq, record.seq) for ref in self.sequences]
            scores_rev = [align_global_score(ref.seq, record.seq.reverse_complement()) for ref in self.sequences]
            if max(scores_fwd) > max(scores_rev):
                idx = scores_fwd.index(max(scores_fwd))
                record.assign_reference(self.sequences[idx].id)

            else:
                idx = scores_rev.index(max(scores_rev))
                record.assign_reference(self.sequences[idx].id)
                record.flag_as_reverse(True)


