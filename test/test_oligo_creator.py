import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import split_guide, oligo_creator

class TestOligoCreator:

    def test_create_oligos(self):
        test_split_guides = [split_guide.SplitGuide(Seq("GGGCGAGGAGCTGTT"), Seq("TGTTCACCG"), Seq("TGTT")),
                             split_guide.SplitGuide(Seq("GCTGCACGCCGTAGG"), Seq("TAGGTCAGGG"), Seq("TAGG")),
                             split_guide.SplitGuide(Seq("GGTGAACCGCATCG"), Seq("ATCGAGCTGA"), Seq("ATCG"))]

        results = oligo_creator.OligoCreator(test_split_guides).create_oligos()


        assert(len(results) == 3)
        assert(str(results[0].name()) == "oligos for [start piece of guide 0] to [end piece of guide 2]")
        assert(str(results[1].name()) == "oligos for [start piece of guide 1] to [end piece of guide 0]")
        assert(str(results[2].name()) == "oligos for [start piece of guide 2] to [end piece of guide 1]")
        assert(str(results[0].forward_oligo) == "forward oligo with overlap ATCG: CACC|GGGCGAGGAGCTGTT|TGAGACCGAGGTCTCA|ATCGAGCTGA")

    def test_create_forward_and_reverse_oligos(self):
        test_split_guide_a = split_guide.SplitGuide(Seq("AAAA"), Seq("GCTGAAAAA"), Seq("GCTG"))
        test_split_guide_b = split_guide.SplitGuide(Seq("GGGG"), Seq("AGGTGGGGG"), Seq("AGGT"))
        results = oligo_creator.OligoCreator([])._create_forward_and_reverse_oligos(test_split_guide_a, test_split_guide_b)

        assert(len(results) == 2)
        assert(str(results[0]) == "forward oligo with overlap AGGT: CACC|AAAA|TGAGACCGAGGTCTCA|AGGTGGGGG")
        assert(str(results[1]) == "reverse oligo with overlap ACCT: AAAC|CCCCCACCT|TGAGACCTCGGTCTCA|TTTT")
        assert(Seq(self._get_raw_sequence(results[0])[4:]).reverse_complement() == Seq(self._get_raw_sequence(results[1])[4:]))

    def _get_raw_sequence(self, oligo):
        return "".join(str(oligo).split(":")[1].strip().split("|"))
