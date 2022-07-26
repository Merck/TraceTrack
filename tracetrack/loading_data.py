import os
import zipfile
from typing import List
from flask import flash
from tracetrack.entities.record import TraceSeqRecord
from tracetrack.entities.errors import ParsingError


def unzip_and_get_sequences(dir: str) -> List[TraceSeqRecord]:
    """
    Process all files in given directory, return list of TraceSeqRecord objects.
    :param dir: str
    :param threshold: int, score threshold for a position to be taken into account
    :param end_threshold: int, score threshold for low-quality end trimming
    :return: list of TraceSeqRecord objects
    """
    filelist = load_all_files(dir)      # list of file names ending with .ab1
    seqlist = []
    for f in filelist:
        try:
            rec = TraceSeqRecord.read(f)
            seqlist.append(rec)
        except ParsingError as error:
            flash(error.message)
    return seqlist


# called by unzip_and_get_sequences
def load_all_files(dir):
    """
    Check all files in directory, unzip zip files, return a list of .ab1 files. Skip invalid files and flash a warning message.
    :param dir: str
    :return: list of file names (str)
    """
    paths = []
    for name in os.listdir(dir):
        if name.lower().endswith('.zip'):
            with zipfile.ZipFile(os.path.join(dir, name), 'r') as zip_ref:
                zip_ref.extractall(dir)
            os.remove(os.path.join(dir, name))
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.lower().endswith('.ab1'):
                paths.append(os.path.join(root, file))
            else:
                flash("File {} has an invalid extension.".format(file))
    # paths.sort()
    return paths




