import pytest

from cargo_oligo_creator import search_offset

class TestSearchOffset:

    def test_get_offsets_for_one(self):
        new_offsets = search_offset.SearchOffset(1)
        all_offsets = new_offsets.get_offsets()
        expected = [[0], [-1], [1], [-2], [2], [-3], [3]]

        assert(all_offsets == expected)

    def test_get_offsets_for_two(self):
        new_offsets = search_offset.SearchOffset(2)
        all_offsets = new_offsets.get_offsets()
        expected = [[0, 0], [0, -1], [0, 1], [0, -2], [0, 2], [0, -3], [0, 3], [-1, 0], [-1, -1], [-1, 1]]

        assert(len(all_offsets) == 49)
        assert(all_offsets[:10] == expected)
        assert(all_offsets[-1] == [3,3])

    def test_get_offsets_for_three(self):
        new_offsets = search_offset.SearchOffset(3)
        all_offsets = new_offsets.get_offsets()
        expected = [[0, 0, 0], [0, 0, -1], [0, 0, 1]]

        assert(len(all_offsets) == 7*7*7)
        assert(all_offsets[:3] == expected)
        assert(all_offsets[-1] == [3,3,3])
