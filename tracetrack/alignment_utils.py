from typing import List, Tuple, Optional
from Bio.Seq import Seq, MutableSeq
from Bio import Align, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import os
import tempfile
import string
from tracetrack.entities.errors import AlignmentError
from tracetrack.entities.record import Record
import parasail
import re

class Settings:
    def __init__(self, threshold=0, end_threshold=0, separate=False):
        self.threshold = threshold
        self.end_threshold = end_threshold
        self.separate = separate


def preserve_lowercase(raw_seq, aligned_seq) -> Tuple[Seq, List[int]]:
    """
    Convert bases in aligned sequence to lowercase or uppercase based on original sequence
    :param raw_seq: original sequence
    :param aligned_seq: aligned sequence with gaps
    :return: aligned sequence with bases converted to lowercase or uppercase based on original sequence
    """
    pos = 0
    mapping = []
    mutable_seq = MutableSeq(aligned_seq)
    for i, base in enumerate(aligned_seq):
        if base == '-':
            continue
        raw_base = raw_seq[pos]
        assert raw_base.upper() == base.upper(), \
            f'Unexpected error, base "{base.upper()}" != "{raw_base.upper()}" on position {pos}:\n' \
            f'    raw: {raw_seq}\n' \
            f'aligned: {aligned_seq}'
        if raw_base.islower():
            mutable_seq[i] = mutable_seq[i].lower()
        else:
            mutable_seq[i] = mutable_seq[i].upper()
        mapping.append(i)
        pos += 1
    return Seq(mutable_seq), mapping


def align_clustalo(ref_record: Record, query_records: List[Record]) -> Tuple[Seq, List[Seq]]:
    """
    Create multiple alignment and return aligned reference and query strings
    :param ref_record: reference record
    :param query_records: query records
    :return: tuple with (aligned reference, list of aligned queries)
    """
    #
    # FIXME set proper gap penalties!
    #
    records = [ref_record] + query_records

    # ensure that all IDs are unique
    ids = [r.id for r in records]
    # TODO
    assert len(set(ids)) == len(ids), f"Found duplicate ID in: {ids}"
    assert ref_record.id, f"Reference ID not found in record: {ref_record}"

    # create tmp input file
    tmp_fd, tmp_path = tempfile.mkstemp()
    result_path = tmp_path + '.aln'
    write_records_to_fasta(records, tmp_fd)

    # generate alignment
    cline = ClustalOmegaCommandline(infile=tmp_path, outfile=result_path, verbose=False, auto=True)
    stdout, stderr = cline()
    if stderr.strip():
        raise AlignmentError(f'Unexpected ClustalO error: {stderr}')

    # parse alignment
    aligned_records = list(AlignIO.read(result_path, "fasta"))
    assert len(aligned_records) == len(records)

    # remove tmp files
    os.remove(tmp_path)
    os.remove(result_path)

    # parse reference and query align strings
    aligned_records_dict = {r.id: r.seq for r in aligned_records}
    ref_aligned_seq = aligned_records_dict[ref_record.id]
    query_aligned_seqs = []
    for r in query_records:
        if r.id not in aligned_records_dict:
            raise AlignmentError(f'Sequence ID "{r.id}" not found in alignment: {aligned_records_dict.keys()}')
        query_aligned_seqs.append(aligned_records_dict[r.id])

    return ref_aligned_seq, query_aligned_seqs


def align_global_score(a, b, gap_open=8, gap_extend=1, matrix=parasail.dnafull):
    """Align two sequences using parwise global alignment with reasonable defaults

    Defaults: BLOSUM62 scoring, gap open = -11, gap extend = -1

    Returns score of alignment
    """
    assert gap_open >= 0
    assert gap_extend >= 0

    result = parasail.nw_trace_striped_16(str(a), str(b), gap_open, gap_extend, matrix)
    return result.score


def translate_codon(seq: str) -> str:
    codon = Seq(seq)
    AA = codon.translate()
    return AA


def write_records_to_fasta(records, filename):
    with open(filename, 'w') as f:
        for record in records:
            f.write('>' + record.id + '\n')
            seq_parts = [str(record.seq[i:min(i+80, len(record.seq))]) for i in range(0, len(record.seq), 80)]
            f.writelines(seq_parts)
            f.write('\n')


def suffix(num):
    """
    Create letter code for numbering insertions. If there is no insertion, return empty string.
    :param num: Number of inserted positions.
    :return: Letter code for insertion or empty string
    """
    res = ""
    if num <= 0:
        return res
    letters = [""] + list(string.ascii_lowercase)
    # use approach similar to number base conversion to generate letter insertion code
    while num > 0:
        modulo = num % 26
        if modulo == 0:
            modulo = 26
            num -= 1
        res = letters[modulo] + res
        num = num // 26
    return res
