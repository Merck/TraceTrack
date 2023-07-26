from typing import List
from collections import Counter


UNKNOWN_INSERTION = '?'


class SequencePosition:
    def __init__(self, ref: str, options: List[str], pos_ref: int, pos_mutated: int, mixed_peak: bool):
        """
        :param ref: reference base
        :param options: list of considered trace bases (no Ns or gaps at edges)
        :param pos_ref: corresponding position in unaligned reference sequence
        :param pos_mutated: corresponding CODING position in mutated result sequence (gaps-deletions are not counted)
        :param mixed_peak: True if there is a mixed peak in the chromatogram at given position
        """
        self.ref = ref
        self.options = options
        self.result, self.confidence, self.coverage = self.vote()
        self.pos_ref = pos_ref
        self.pos_mutated = pos_mutated
        self.mixed_peak = mixed_peak

    def __len__(self):
        return len(self.options)

    def is_mutated(self):
        return self.ref.upper() != self.result.upper()

    def has_coverage(self):
        return len(self.options) > 0

    def is_indel(self):
        return self.ref == "-" or self.result == "-"

    def vote(self):
        """
        Choose the correct base for given SequencePosition. If no options are available, return reference.
        If all options contain the same mutation, return that. If at least one option matches reference, return that option.
        If the options are different and do not match reference, return reference.
        Coverage is given by number of options that match the returned base. Confidence is coverage / number of options.
        :return: reference base (string), confidence (float), coverge (int)
        """
        num_of_options = len(self.options)
        if num_of_options == 0:
            return self.ref, 0, 0  # there are no options -> return reference
        occurrence_count = Counter(self.options)
        most_common_letter, most_common_count = occurrence_count.most_common(1)[0]
        if most_common_count == num_of_options:  # all reads contain the same letter -> return that
            result = most_common_letter.lower() if self.ref.islower() else most_common_letter.upper()
            return result, 1, num_of_options
        # return reference coverage reduced to what matches
        coverage = sum(o.upper() == self.ref.upper() for o in self.options)  # number of reads that match reference
        # what should the confidence be here? Include all letters, or just the ones used for coverage?
        confidence = coverage / num_of_options
        # if the mutation is an insert (gap in reference) and we don't agree, return '?'
        if self.ref == '-':
            return UNKNOWN_INSERTION, 0, 0
        return self.ref, confidence, coverage

    def get_codon_num(self, mutated: bool) -> int:
        """Get the position number of the codon at given position in the sequence."""
        if mutated:
            return self.pos_mutated // 3
        else:
            return self.pos_ref // 3

    def mixed_peak_css_class(self):
        """Return string for formatting html sequence display."""
        if self.mixed_peak:
            return "mixed_peak"
        else:
            return ""

