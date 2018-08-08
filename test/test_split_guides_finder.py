import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import guide
from cargo_oligo_creator import split_guides_finder

class TestSplitGuidesFinder:

    def test_find_split_guides(self):
        sequences = ["GGGCGAGGAGCTGTTCACCG", "GCTGCACGCCGTAGGTCAGGG", "GGTGAACCGCATCGAGCTGA"]
        test_guides = [guide.Guide(x) for x in sequences]
        results = split_guides_finder.SplitGuidesFinder(test_guides).find_split_guides()

        assert(len(results) == 3)
        # all guides should've split at the halfway point
        first_split_guide = [results[0].first_part, results[0].second_part, results[0].overlap]
        expected_first_split_guide = [Seq("GGGCGAGGAGCTGTT"), Seq("TGTTCACCG"), Seq("TGTT")]
        assert(first_split_guide == expected_first_split_guide)
        second_split_guide = [results[1].first_part, results[1].second_part, results[1].overlap]
        expected_second_split_guide = [Seq("GCTGCACGCCGTAGG"), Seq("TAGGTCAGGG"), Seq("TAGG")]
        assert(second_split_guide == expected_second_split_guide)
        third_split_guide = [results[2].first_part, results[2].second_part, results[2].overlap]
        expected_third_split_guide = [Seq("GGTGAACCGCATCG"), Seq("ATCGAGCTGA"), Seq("ATCG")]
        assert(third_split_guide == expected_third_split_guide)
