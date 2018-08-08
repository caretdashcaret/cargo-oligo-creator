from .split_guide import SplitGuide

class GuideSplitter:
    OVERLAP = 4

    def __init__(self, guides):
        self.guides = guides

    def split(self, split_configuration):
        splits = []
        for index, guide in enumerate(self.guides):
            split_index = split_configuration.get_split_index_for_guide(index)
            split_guide = SplitGuide(
                first_part=self._split_first_part(guide, split_index),
                second_part=self._split_second_part(guide, split_index),
                overlap=self._find_overlap(guide, split_index))
            splits.append(split_guide)
        return splits

    def _split_first_part(self, guide, split_index):
        return guide.sequence[:(split_index + 1 + self.OVERLAP)]

    def _split_second_part(self, guide, split_index):
        return guide.sequence[split_index + 1:]

    def _find_overlap(self, guide, split_index):
        return guide.sequence[split_index + 1: (split_index + 1 + self.OVERLAP)]
