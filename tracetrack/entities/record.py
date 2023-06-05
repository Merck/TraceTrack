from Bio import SeqRecord, SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import FeatureLocation
import os
import numpy as np
from tracetrack.entities.errors import ReferenceError


TRACE_CHANNELS = {
    'G': 'DATA9',
    'A': 'DATA10',
    'T': 'DATA11',
    'C': 'DATA12'
}


class Feature:
    def __init__(self, type, location, label, fid):
        assert isinstance(location.start, int), 'Location start should be an integer'
        assert isinstance(location.end, int), 'Location end should be an integer'
        self.type = type
        self.label = label
        self.location = location
        self.id = fid

    def get_format(self):
        return "default-feature"

    def is_coding(self):
        return self.type == "CDS"


class CDS(Feature):
    def __init__(self, location, label, fid):
        assert isinstance(location.start, int), 'Location start should be an integer'
        assert isinstance(location.end, int), 'Location end should be an integer'
        super().__init__("CDS", location, label, fid)

    def get_format(self):
        return "cds"

    def is_coding(self):
        return True


class Record:
    def __init__(self, sequence: Seq, id: str, features=[]):
        self.seq = sequence
        self.id = id
        self.features = features

    @classmethod
    def read_from_genbank(cls, seq_record: SeqRecord):
        features = []
        cds = False
        i = 0
        for feature in seq_record.features:
            try:
                label = feature.qualifiers['label'][0]
            except KeyError as e:
                label = "none"
            if feature.type == "CDS":
                features.append(CDS(feature.location, label, i))
                cds = True
            else:
                features.append(Feature(feature.type, feature.location, label, i))
            i += 1
        if not cds:
            raise ReferenceError("Each reference sequence must contain at least one CDS feature.")
        return cls(seq_record.seq, seq_record.id, features)

    @classmethod
    def read_from_fasta(cls, seq_record: SeqRecord):
        loc = FeatureLocation(0, len(seq_record.seq))
        return cls(seq_record.seq, seq_record.id, [CDS(loc, "none", 0)])

    @classmethod
    def read_from_csv(cls, row):
        loc = FeatureLocation(0, len(row[1]))
        return cls(Seq(row[1]), row[0], [CDS(loc, "none", 0)])

    def __len__(self):
        return len(self.seq)


class TraceSeqRecord(Record):
    def __init__(self, seq, quality, mixed_fraction, id=None, traces=None, base_locations=None, reference=None,
                 reverse: bool = None):
        assert not isinstance(seq, str), 'Sequence should be a Seq object'
        super().__init__(seq, id)
        self.traces = traces
        self.base_locations = base_locations
        self.quality = quality
        self.mixed_peaks = self.find_mixed_peaks(mixed_fraction)
        # should the original sequence be stored because of finding mixed peaks? Or just ignore al N's...
        self.reference = reference
        self.reverse = reverse

    @classmethod
    def read(cls, path, mixed_fraction):
        """
        Take an .ab1 file, return TraceSeqRecord object. Raise exception when the file name format is wrong.
        :param path: str
        :param threshold: int, score threshold for a position to be taken into account
        :param end_threshold: int, score threshold for low-quality end trimming
        :return: TraceSeqRecord
        """
        name_w_ext = os.path.basename(path)  # removes path
        name = os.path.splitext(name_w_ext)[0]

        # load record with sequence
        record = SeqIO.read(path, 'abi')

        # get trace values for all bases
        traces = {base: list(record.annotations['abif_raw'][channel]) for base, channel in TRACE_CHANNELS.items()}

        # get locations in trace array for all bases
        base_locations = list(record.annotations['abif_raw']["PLOC1"])

        return cls(record.seq, record.letter_annotations['phred_quality'], mixed_fraction=mixed_fraction, id=name,
                   traces=traces, base_locations=base_locations)

    def filter_sequence_by_quality(self, threshold, end_threshold):
        """
        replace low-quality ends and all bases under quality threshold with Ns, return copy of object
        """
        beg, end = find_quality_ends(self.quality, end_threshold)
        mutable_seq = MutableSeq(self.seq)
        for pos, q in enumerate(self.quality):
            if q < threshold or pos < beg or pos >= end:
                mutable_seq[pos] = 'N'
        return TraceSeqRecord(
            Seq(mutable_seq),
            self.quality,
            id=self.id,
            traces=self.traces,
            base_locations=self.base_locations,
            reference=self.reference,
            reverse=self.reverse
        )

    def reverse_complement(self, **kwargs):
        num_locations = len(self.traces['A'])
        return TraceSeqRecord(
            self.seq.reverse_complement(**kwargs),
            self.quality,
            id=self.id,
            traces={str(Seq(base).reverse_complement()): values[::-1] for base, values in self.traces.items()},
            base_locations=[num_locations - i - 1 for i in self.base_locations[::-1]],
            reverse=False,
            reference=self.reference
        )

    def has_base_above_threshold(self):
        for i in range(len(self.seq)):
            if self.seq[i] != "N":
                return True
        return False

    def quality_counts(self):
        counts = np.zeros(100)
        for quality in self.quality:
            counts[quality] += 1
        return counts

    def assign_reference(self, refid):
        self.reference = refid

    def flag_as_reverse(self, reverse):
        self.reverse = reverse

    def area_under_peak(self, start: int, end: int, base: str):
        """
        Return area under the curve in a given trace between indices start and end (including both)
        :param start: start index for curve region
        :param end: end index
        :param base: specifies trace (A, C, T or G)
        :return: number corresponding to (approximate) area under curve
        """
        area = sum(self.traces[base.upper()][start:end+1])
        return area

    def peak_borders(self, i: int):
        """
        Determine positions in the trace that correspond to a single chromatogram peak
        :param i: position in sequence
        :return: border positions in trace
        """
        pos = self.base_locations[i]
        if i == 0:
            prev = -pos
        else:
            prev = self.base_locations[i - 1]
        if i == len(self.base_locations) - 1:
            end = len(self.traces["A"])
            start = (pos - prev) // 2
        else:
            next = self.base_locations[i + 1]
            width = (next - prev) // 4
            start = max(0, pos - width + 2)
            end = min(pos + width - 2, len(self.traces['A']))
        return start, end

    def find_mixed_peaks(self, fraction: float = 0.15):
        """
        For each position of the sequence, determine if the peak in the chromatogram is "mixed".
        Disregard mixed signals in regions with low signal to noise ratio (generally bad quality region)
        :param fraction: Threshold for ratio of secondary to primary peak area that is considered mixed
        :return: List of sequence positions with mixed signal
        """
        # precalculate signal to noise ratio per position
        stn = [self.signal_to_noise(i) for i in range(len(self.base_locations))]
        avg_stn = sum(stn) / len(stn)
        mixed_peaks = []

        for i, pos in enumerate(self.base_locations):
            stn_local = stn[i-10:i] + stn[i+1:i+10]
            signal_to_noise = sum(stn_local) / 20

            threshold = max(25, avg_stn * 1.35)
            if signal_to_noise < threshold:
                continue
                # bad StN ratio -> disregard potential mixed positions

            base = self.seq[i]
            start, end = self.peak_borders(i)
            areas = {base: self.area_under_peak(start, end, base) for base in self.traces.keys()}
            peaks = {base: values[pos] for base, values in self.traces.items()}
            if base != "N":
                main_peak = peaks[base.upper()]
                for letter, area in areas.items():
                    # check for both area and height of peak
                    if base != letter and area > (areas[base.upper()] * fraction) and peaks[letter] > (main_peak * fraction) \
                            and self.is_concave(pos, letter):
                        mixed_peaks.append(i)
        return mixed_peaks

    def signal_to_noise(self, i: int):
        """
        Calculate the signal to noise ratio of a position as the ratio of the primary peak area to the sum of all other
        signal areas
        :param i: position in sequence
        :return: Signal to noise ratio
        """
        base = self.seq[i]
        if base.upper() == 'N':
            return 1
        start, end = self.peak_borders(i)
        #base = self.seq[i]
        areas = {base: self.area_under_peak(start, end, base) for base in self.traces.keys()}
        primary = areas[base.upper()]
        secondary = sum([areas[letter.upper()] for letter in areas.keys() if letter != base])
        if secondary == 0:
            # a sufficiently high number for positions with no signal from secondary bases
            return 35
        return primary / secondary

    def check_area_quality(self, i: int, window: int = 10):
        ratio = 0
        for base in range(i - window, i + window):
            if base != i:
                ratio += self.signal_to_noise(base)
        return ratio / window * 2

    def is_concave(self, pos: int, base: str):
        """
        Determine whether the trace at given position is concave
        :param pos: position in the trace (not in the sequence)
        :param base: base - which trace are we looking at
        :return: bool (true if concave)
        """
        x1 = pos - 2
        x2 = pos
        x3 = pos + 2
        func = self.traces[base.upper()]
        left = (func[x2] - func[x1]) / (x2 - x1)
        right = (func[x3] - func[x1]) / (x3 - x1)
        return left >= right


def find_quality_ends(scores, threshold):
    """
    Search from both ends of the sequence until three consecutive bases with quality exceeding threshold are found,
    return positions where this occurs
    :param scores: list of per base qualities
    :param threshold: end trimming quality threshold
    :return: indices, where sequence is trimmed
    """
    i = -1
    beg = i
    found = False
    high_scores = 0
    while i < len(scores) - 1 and not found:
        i += 1
        if scores[i] >= threshold:
            high_scores += 1
        else:
            high_scores = 0
            beg = i
        if high_scores >= 3:
            found = True
    i = len(scores)
    end = len(scores)
    found = False
    high_scores = 0
    while not found and i > 0:
        i -= 1
        if scores[i] >= threshold:
            high_scores += 1
        else:
            high_scores = 0
            end = i
        if high_scores >= 3:
            found = True
    return beg + 1, end
