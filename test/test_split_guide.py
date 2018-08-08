import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import split_guide

class TestSplitGuide:

    def test_create_split_guide(self):
        first_part = Seq("AAATT")
        second_part = Seq("TTGCG")
        overlap = Seq("TT")
        new_split_guide = split_guide.SplitGuide(
            first_part=first_part,
            second_part=second_part,
            overlap=overlap)

        assert(new_split_guide.first_part == first_part)
        assert(new_split_guide.second_part == second_part)
        assert(new_split_guide.overlap == overlap)
