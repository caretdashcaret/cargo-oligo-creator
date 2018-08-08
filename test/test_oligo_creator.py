import pytest

from cargo_oligo_creator import guide, oligo_creator

class TestOligoCreator:

    def test_create_oligos(self):
        sequences = ["GGGCGAGGAGCTGTTCACCG", "GCTGCACGCCGTAGGTCAGGG", "GGTGAACCGCATCGAGCTGA"]
        test_guides = [guide.Guide(x) for x in sequences]
        results = oligo_creator.OligoCreator(test_guides).create_oligos()

        assert(len(results) == 4)
        assert(str(results[0]) == "Forward oligo for [end of guide 0] to [start of guide 1] with overlap TAGG: CACC|TGTTCACCG|TGAGACCGAGGTCTCA|GCTGCACGCCGTAGG")
        assert(str(results[1]) == "Reverse oligo for [end of guide 0] to [start of guide 1] with overlap ACCG: AAAC|CCTACGGCGTGCAGC|ACTCTGGCTCCAGAGT|CGGTGAACA")
        assert(str(results[2]) == "Forward oligo for [end of guide 1] to [start of guide 2] with overlap ATCG: CACC|TAGGTCAGGG|TGAGACCGAGGTCTCA|GGTGAACCGCATCG")
        assert(str(results[3]) == "Reverse oligo for [end of guide 1] to [start of guide 2] with overlap AGGG: AAAC|CGATGCGGTTCACC|ACTCTGGCTCCAGAGT|CCCTGACCTA")


    def test_create_oligos(self):
        sequences = ["GGGCGAGGAGCTGTTCACCG", "GCTGCACGCCGTAGGTCAGGG", "GGTGAACCGCATCGAGCTGA", "GGTGTTCTGCTGGTAGTGGT"]
        test_guides = [guide.Guide(x) for x in sequences]
        results = oligo_creator.OligoCreator(test_guides).create_oligos()

        assert(len(results) == 6)
