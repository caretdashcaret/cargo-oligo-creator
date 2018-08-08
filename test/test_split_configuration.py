import pytest

from cargo_oligo_creator import split_configuration

class TestSplitConfiguration:

    def test_get_split_index_for_guide(self):
        new_config = split_configuration.SplitConfiguration([20])

        assert(new_config.get_split_index_for_guide(0) == 20)
