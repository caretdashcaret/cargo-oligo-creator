import pytest

from cargo_oligo_creator import oligos_display

class TestOligosDisplay:

    def test_oligos_display(self):
        test_display = oligos_display.OligosDisplay(
            forward_oligo=1,
            reverse_oligo=2,
            pre_construct_guide_index = 2,
            post_construct_guide_index = 10)

        assert(test_display.name() == "oligos for [start piece of guide 2] to [end piece of guide 10]")
        assert(test_display.forward_oligo == 1)
        assert(test_display.reverse_oligo == 2)
