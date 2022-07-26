import pytest
from tracetrack.entities.alignment import Alignment
from Bio.Seq import Seq
from tracetrack.entities.record import TraceSeqRecord, Record, CDS, FeatureLocation
from tracetrack.entities.sequence_position import SequencePosition
from tracetrack.alignment_utils import align_global_score
from tracetrack.entities.reference_db import ReferenceDb
from tempfile import TemporaryDirectory
import shutil
from tracetrack.server_utils import unzip_and_get_sequences
from tracetrack.tasks import align_sequences_celery, AlignmentTaskResult
from tracetrack.server_utils import schedule_tasks, scheduler
from tracetrack.entities.scheduler import use_scheduler
from time import sleep
from typing import List
use_scheduler('simple')


class TestAlignmentCase:
    def test_vote(self):
        seq_pos = SequencePosition("T", ["A", "A", "G"], 34, 33, False)
        assert seq_pos.vote() == ("T", 0, 0)  # result, confidence, coverage
        seq_pos = SequencePosition("T", ["T", "G", "A"], 34, 34, False)
        assert seq_pos.vote() == ("T", 1 / 3, 1)
        seq_pos = SequencePosition("T", ["C", "C", "C"], 34, 34, False)
        assert seq_pos.vote() == ("C", 1, 3)
        seq_pos = SequencePosition("T", [], 34, 34, False)
        assert seq_pos.vote() == ("T", 0, 0)

    """
          0123456789
    ref:  CCAATGATGCGG
    seq1: NNAATNTTGCGNN
    seq2: NNAATGTAGCGNN
    mut:  CCAATGTTGCGG
    """

    def test_alignment(self, mocker):
        mocker.patch('tracetrack.entities.record.TraceSeqRecord.find_mixed_peaks', return_value=[])
        ref = "CCAATGATGCGG"
        ref_rec = Record(Seq(ref), "1", [CDS(FeatureLocation(0, 12, 1), "cds", 0)])
        seq_rec1 = TraceSeqRecord(Seq("NNAATNTTGCGNN"), [], "1-tracefile1.ab1")
        seq_rec2 = TraceSeqRecord(Seq("NNAATGTAGCGNN"), [], "1-tracefile2.ab1")
        alignment1 = Alignment.align_multiple(ref_rec, [seq_rec1, seq_rec2], "pop1")
        mutated = [pos.result for pos in alignment1.positions]

        # num	ref	s1	s2	mut
        #	0	C	N	N	C
        #	1	C	N	N	C
        #	2	A	A	A	A
        #	3	A	A	A	A
        #	4	T	T	T	T
        #	5	G	N	G	G => keeping original
        #	6	A	T	T	T => missense mutation
        #	7	T	T	A	T => ???? mutation
        #	8	G	G	G	G
        #	9	C	C	C	C
        #		G	G	G	G
        #		G	N	N	G
        #			N	N

        assert mutated == list('CCAATGTTGCGG')
        assert alignment1.mutations == 1

        seqlist = alignment1.positions
        assert seqlist[3].ref == "A"
        assert seqlist[4].options == ["T", "T"]

        assert alignment1.get_codon_reference(6) == "ATG"
        assert alignment1.get_codon_mutated(7) == "TTG"

        assert alignment1.positions[0].coverage == 0
        assert alignment1.positions[3].coverage == 2

        assert alignment1.get_base_at_position(3, True) == "A"
        assert alignment1.get_base_at_position(3, False) == "A"
        assert alignment1.get_base_at_position(6, True) == "T"
        assert alignment1.get_base_at_position(6, False) == "A"

        assert alignment1.is_same_aminoacid(3)
        assert not alignment1.is_same_aminoacid(6)
        assert not alignment1.is_same_aminoacid(7)

        assert alignment1.mutation_type(6) == 'mis_mut'
        assert alignment1.mutation_type(3) == 'no_mut'

    def test_alignment_with_indels(self, mocker):
        mocker.patch('tracetrack.entities.record.TraceSeqRecord.find_mixed_peaks', return_value=[])
        ref = "atggagtggagctgggtgttcctgttcttcctgagcgtgaccaccggcgtgcacagcCAAGTGCAGCTGCAAGAGAGCGGACC"
        ref_rec = Record(Seq(ref), "1", [CDS(FeatureLocation(0, 90, 1), "cds", 0)])
        seq_rec1 = TraceSeqRecord(Seq("CACCATGGAGTGGAGCTGGGTGTTCCTTTCTTCCTGAGCGTGACCACCGGCGTGCACAGCCAAAGTGCAGCTGCAAGAGAGCGGACCCGG"), [], "1-tracefile1.ab1")
        seq_rec2 = TraceSeqRecord(Seq("CACCATGGAGTGGAGCTGGGTGTTCCTTTCTTCCTGAGCGTGACCACCGGCGTGCACAGCCAAGGTGCCGCTGCAAGAGAGCGGACCCGG"), [], "1-tracefile2.ab1")
        alignment = Alignment.align_multiple(ref_rec, [seq_rec1, seq_rec2], "pop1")
        mutated = [pos.result for pos in alignment.positions]

        """
        What I'd expect:
        ref: atggagtggagctgggtgttcctgttcttcctgagcgtgaccaccggcgtgcacagcCAA-GTGCAGCTGCAAGAGAGCGGACC
        tr1: ATGGAGTGGAGCTGGGTGTTCCT-TTCTTCCTGAGCGTGACCACCGGCGTGCACAGCCAAAGTGCAGCTGCAAGAGAGCGGACC
        tr2: ATGGAGTGGAGCTGGGTGTTCCT-TTCTTCCTGAGCGTGACCACCGGCGTGCACAGCCAAGGTGCCGCTGCAAGAGAGCGGACC
        mut: atggagtggagctgggtgttcct-ttcttcctgagcgtgaccaccggcgtgcacagcCAA?GTGCAGCTGCAAGAGAGCGGACC
        
        What Clustal actually does:
        ref: atggagtggagctgggtgttcctgttcttcctgagcgtgaccaccggcgtgcacagcCAA-GTGCAGCTGCAAGAGAGCGGACC
        tr1: ATGGAGTGGAGCTGGGTGTT-CCTTTCTTCCTGAGCGTGACCACCGGCGTGCACAGCCAAAGTGCAGCTGCAAGAGAGCGGACC
        tr2: ATGGAGTGGAGCTGGGTGTT-CCTTTCTTCCTGAGCGTGACCACCGGCGTGCACAGCCAAGGTGCCGCTGCAAGAGAGCGGACC
        mut: atggagtggagctgggtgtt-cctttcttcctgagcgtgaccaccggcgtgcacagcCAA?GTGCAGCTGCAAGAGAGCGGACC
        """
        assert mutated == list('atggagtggagctgggtgtt-cctttcttcctgagcgtgaccaccggcgtgcacagcCAA?GTGCAGCTGCAAGAGAGCGGACC')

        assert alignment.get_base_at_position(3, True).upper() == "G"
        assert alignment.get_base_at_position(3, False).upper() == "G"
        assert alignment.get_base_at_position(20, True) == "-"

        assert alignment.mutation_type(6) == 'no_mut'
        assert alignment.mutation_type(22) == 'frameshift'

        assert alignment.frameshifts == 64
        assert alignment.get_base_at_position(63, True) == "G"

    def test_align(self, mocker):
        mocker.patch('tracetrack.entities.record.TraceSeqRecord.find_mixed_peaks', return_value=[])
        ref = "CCAATGATGCGG"
        ref_rec = Record(Seq(ref), "1", [CDS(FeatureLocation(0, 90, 1), "cds", 0)])
        seq_rec1 = TraceSeqRecord(Seq("NNAATNTTGCGNN"), [], "1-tracefile1.ab1")
        alignment = Alignment.align_multiple(ref_rec, [seq_rec1], "pop1")

        assert str(alignment.aligned_traces[0].aligned_seq) == "NNAATNTTGCGN"
        assert not alignment.stop_count
        assert not alignment.kozak_count

    @pytest.mark.parametrize("reference, trace, target_score",
                             [
                                 (
                                         "AAAACGCTACGC",  # reference
                                         "AAACGCTCGC",  # trace
                                         34
                                 ),
                                 (
                                         "AAAACGCTACGC",  # reference
                                         "TTAAAACGCTACGGG",  # trace
                                         34
                                 ),
                                 (
                                         "AAAACGCTACGC",  # reference
                                         "TTAACGCTACG",  # trace
                                         29
                                 )
                             ])
    def test_parasail_aligner(self, reference, trace, target_score):
        score = align_global_score(reference, trace)
        assert score == target_score

    def test_alignment_complete(self, mocker):
        mocker.patch('tracetrack.server_utils.flash', return_value="flash")
        mocker.patch('tracetrack.entities.reference_db.flash', return_vaue="flash_db")
        with TemporaryDirectory() as tmp:
            shutil.copy2('./resources/test_trace_files/crossarchus_2_missense.ab1', tmp)
            shutil.copy2('./resources/test_trace_files/urva_javanica_4_deletion.ab1', tmp)
            seqlist = unzip_and_get_sequences(tmp)
        db = ReferenceDb.read_file("./resources/example_sheet.xlsx")
        for record in seqlist:
            db.match_trace_to_ref(record)
        task_id = schedule_tasks([seqlist], ["pop1"], db, False, 45, 0)
        while not scheduler.are_results_ready(task_id):
            sleep(2)
        results: List[AlignmentTaskResult] = scheduler.get_results(task_id, timeout=10)
        errors = [str(r) for r in results if isinstance(r, Exception)]
        results = [r for r in results if isinstance(r, AlignmentTaskResult)]

        alignments = [res.alignment for res in results]
        warnings = errors + [warn for res in results for warn in res.warnings]

        assert len(alignments) == 2
        alignments.sort(key=lambda x: x.aligned_reference.record.id)
        alignment = alignments[0]
        assert alignment.aligned_reference.record.id == "2_missense"
        assert alignment.mutations == 1
        assert len(alignment.aligned_traces) == 1
        assert alignment.mis_mutations == 1
        assert alignment.positions[4].coverage == 0

        alignment = alignments[1]
        assert alignment.aligned_reference.record.id == "4_deletion"
        assert len(alignment.aligned_traces) == 1
        assert alignment.mis_mutations == 8
        assert alignment.nons_mutations == 0
        assert alignment.positions[135].result == "a"



