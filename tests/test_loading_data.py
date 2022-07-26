from tracetrack.entities.record import TraceSeqRecord
from tracetrack.entities.reference_db import ReferenceDb
from tracetrack.server_utils import unzip_and_get_sequences, load_all_files
from tempfile import TemporaryDirectory
from distutils.dir_util import copy_tree
import shutil
from Bio.Seq import Seq


class TestLoadingDataCase:
    def test_db(self, mocker):
        mocker.patch('tracetrack.entities.reference_db.flash', return_vaue="flash_db")
        filename = "./resources/test_references.xlsx"
        db = ReferenceDb.read_file(filename)
        assert db.get_ref_sequence("3").seq == "GCGTAAGTTCGC"

    def test_align_sequences(self, mocker):
        mocker.patch('tracetrack.entities.record.TraceSeqRecord.find_mixed_peaks', return_value=[])
        mocker.patch('tracetrack.entities.reference_db.flash', return_vaue="flash_db")
        seq_rec1 = TraceSeqRecord(Seq("GAGCCCAGATTACCAGAGCGTGACACCTAAGGAGAAGTCACTATTCACATGTAGGGCTAGGGCTTA"), [], "4-tracefile_a")
        seq_rec2 = TraceSeqRecord(Seq("GAGCCCAGATTACCAGAGCGTGACACCTAAGGAGAAGTCACTATTCACATGTAGGGCTAGGGCTTA"), [], "4-tracefile_b")
        #                         ref: GAGCCCAGATTTCCAGAGCGTGACACCTAAGGAGAAGTCACTATCACATGTAGGGCTAGGGCTTA
        sequences = [seq_rec1, seq_rec2]
        filename = "./../resources/test_references.xlsx"
        db = ReferenceDb.read_file(filename)
        for record in sequences:
            db.match_trace_to_ref(record)
        alignment, warnings = db.align_sequences(sequences, seq_rec1.reference, "pop1")

        assert alignment.mutations == 23
        assert alignment.mis_mutations == 7
        assert alignment.frameshifts == 14
        assert alignment.silent_mutations == 1

    def test_load_all_files(self, mocker):
        mocker.patch('tracetrack.server_utils.flash', return_value="flash")
        directory = "./resources/test_trace_files"
        with TemporaryDirectory() as tmp:
            copy_tree(directory, tmp)
            file_list = load_all_files(tmp)
            file_list.sort()
        assert file_list == [tmp+"/crossarchus_2_missense.ab1", tmp+"/test_file1.ab1", tmp+"/test_file2.ab1", tmp+"/urva_javanica_4_deletion.ab1"]


def test_loading_data(mocker):
    mocker.patch('tracetrack.server_utils.flash', return_value="flash")
    mocker.patch('tracetrack.entities.reference_db.flash', return_vaue="flash_db")
    with TemporaryDirectory() as tmp:
        shutil.copy2('../resources/test_trace_files/Archive.zip', tmp)
        shutil.copy2('../resources/test_trace_files/crossarchus_2_missense.ab1', tmp)
        shutil.copy2('../resources/test_trace_files/urva_javanica_4_deletion.ab1', tmp)
        seqlist = unzip_and_get_sequences(tmp)
        db = ReferenceDb.read_file("../resources/test_references.xlsx")
        for record in seqlist:
            db.match_trace_to_ref(record)
        alignment1, warnings1 = db.align_sequences([seqlist[0]], seqlist[0].reference, "pop1")
        alignment2, warnings2 = db.align_sequences([seqlist[1]], seqlist[1].reference, "pop1")

        assert alignment1 is not None
        assert alignment2 is not None



