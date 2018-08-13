import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import split_guide, oligo_creator

class TestOligoCreator:

    # def test_create_oligos(self):
    #     test_split_guides = [split_guide.SplitGuide("TGTTCACCG", "GCTGCACGCCGTAGG", "GCTG")]
    #     results = oligo_creator.OligoCreator(test_split_guides).create_oligos()

    # #     class SplitGuide:

    # # def __init__(self, first_part, second_part, overlap):
    # #     self.first_part = first_part
    # #     self.second_part = second_part
    # #     self.overlap = overlap


    #     assert(len(results) == 2)
    #     assert(str(results[0]) == "Forward oligo for [end of guide 0] to [start of guide 1] with overlap TAGG: CACC|TGTTCACCG|TGAGACCGAGGTCTCA|GCTGCACGCCGTAGG")
    #     assert(str(results[1]) == "Reverse oligo for [end of guide 0] to [start of guide 1] with overlap ACCG: AAAC|CCTACGGCGTGCAGC|ACTCTGGCTCCAGAGT|CGGTGAACA")

    def test_create_forward_and_reverse_oligos(self):
        test_split_guide_a = split_guide.SplitGuide(Seq("AAAA"), Seq("GCTGAAAAA"), Seq("GCTG"))
        test_split_guide_b = split_guide.SplitGuide(Seq("GGGG"), Seq("AGGTGGGGG"), Seq("AGGT"))
        results = oligo_creator.OligoCreator([])._create_forward_and_reverse_oligos(test_split_guide_a, test_split_guide_b)

        assert(len(results) == 2)
        assert(str(results[0]) == "Forward oligo with overlap AGGT: CACC|AAAA|TGAGACCGAGGTCTCA|AGGTGGGGG")
        assert(str(results[1]) == "Reverse oligo with overlap ACCT: AAAC|CCCCCACCT|TGAGACCTCGGTCTCA|TTTT")
        assert(Seq(self._get_raw_sequence(results[0])[4:]).reverse_complement() == Seq(self._get_raw_sequence(results[1])[4:]))

    def _get_raw_sequence(self, oligo):
        return "".join(str(oligo).split(":")[1].strip().split("|"))
