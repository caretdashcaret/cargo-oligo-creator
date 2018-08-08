from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class Guide:

    def __init__(self, sequence_string):
        self.sequence = Seq(sequence_string.upper(), generic_dna)
        self.center = self._find_center()

    def _find_center(self):
        return int(self._find_length()/2)

    def _find_length(self):
        return len(self.sequence)
