import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import guide

class TestGuide:

    def test_create_guide(self):
        sequence = "AATTaattCGCGCG"
        test_guide = guide.Guide(sequence)

        assert(test_guide.sequence == Seq("AATTAATTCGCGCG"))
        assert(test_guide.center == 7)
