import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import oligo

class TestOligo:

    def test_create_forward_oligo(self):
        test_oligo = oligo.Oligo(
            is_forward = True,
            overlap = Seq("AAAA"),
            guide_end = 5,
            guide_start = 4,
            start = Seq("TATC"),
            pre_construct_piece = Seq("GGGG"),
            construct = Seq("CCCC"),
            post_construct_piece = Seq("GGGGG")
        )

        assert(test_oligo.formatted_sequence == "TATC|GGGG|CCCC|GGGGG")
        assert(str(test_oligo) == "Forward oligo for [end of guide 5] to [start of guide 4] with overlap AAAA: TATC|GGGG|CCCC|GGGGG")

    def test_create_reverse_oligo(self):
        test_oligo = oligo.Oligo(
            is_forward = False,
            overlap = Seq("AAAA"),
            guide_end = 2,
            guide_start = 3,
            start = Seq("TATC"),
            pre_construct_piece = Seq("GGGG"),
            construct = Seq("CCCC"),
            post_construct_piece = Seq("GGGGG")
        )

        assert(test_oligo.formatted_sequence == "TATC|GGGG|CCCC|GGGGG")
        assert(str(test_oligo) == "Reverse oligo for [end of guide 2] to [start of guide 3] with overlap AAAA: TATC|GGGG|CCCC|GGGGG")
