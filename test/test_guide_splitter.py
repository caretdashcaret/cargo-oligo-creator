import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import guide, guide_splitter, split_configuration

class TestGuideSplitter:

    def test_split(self):
        guides = [guide.Guide("AAAATTCCCCGG")]
        config = split_configuration.SplitConfiguration([3])
        splitter = guide_splitter.GuideSplitter(guides)
        splits = splitter.split(config)

        assert(len(splits) == 1)
        assert(splits[0].first_part == Seq("AAAATTCC"))
        assert(splits[0].second_part == Seq("TTCCCCGG"))
        assert(splits[0].overlap == Seq("TTCC"))
