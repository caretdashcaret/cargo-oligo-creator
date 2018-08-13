import pytest

from Bio.Seq import Seq
from cargo_oligo_creator import oligo

class TestOligo:

    def test_create_forward_oligo(self):
        test_oligo = oligo.Oligo(
            is_forward = True,
            overlap = Seq("AAAA"),
            start = Seq("TATC"),
            pre_construct_piece = Seq("GGGG"),
            construct = Seq("CCCC"),
            post_construct_piece = Seq("GGGGG")
        )

        assert(test_oligo.raw_sequence() == "TATCGGGGCCCCGGGGG")
        assert(test_oligo.formatted_sequence() == "TATC|GGGG|CCCC|GGGGG")
        assert(test_oligo.display_formatted() == "forward oligo with overlap AAAA: TATC|GGGG|CCCC|GGGGG")
        assert(test_oligo.display_raw() == "forward oligo with overlap AAAA: TATCGGGGCCCCGGGGG")
        assert(test_oligo.direction() == "forward")
        assert(str(test_oligo) == "forward oligo with overlap AAAA: TATC|GGGG|CCCC|GGGGG")

    def test_create_reverse_oligo(self):
        test_oligo = oligo.Oligo(
            is_forward = False,
            overlap = Seq("AAAA"),
            start = Seq("TATC"),
            pre_construct_piece = Seq("GGGG"),
            construct = Seq("CCCC"),
            post_construct_piece = Seq("GGGGG")
        )

        assert(test_oligo.raw_sequence() == "TATCGGGGCCCCGGGGG")
        assert(test_oligo.formatted_sequence() == "TATC|GGGG|CCCC|GGGGG")
        assert(test_oligo.display_formatted() == "reverse oligo with overlap AAAA: TATC|GGGG|CCCC|GGGGG")
        assert(test_oligo.display_raw() == "reverse oligo with overlap AAAA: TATCGGGGCCCCGGGGG")
        assert(test_oligo.direction() == "reverse")
        assert(str(test_oligo) == "reverse oligo with overlap AAAA: TATC|GGGG|CCCC|GGGGG")
