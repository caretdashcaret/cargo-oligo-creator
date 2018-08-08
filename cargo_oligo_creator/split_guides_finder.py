from .guide_splitter import GuideSplitter
from .search_offset import SearchOffset
from .split_configuration import SplitConfiguration

class SplitGuidesFinder:

    def __init__(self, guides):
        self.guides = guides

    def find_split_guides(self):
        guide_splitter = GuideSplitter(self.guides)
        offsets = SearchOffset(len(self.guides)).get_offsets()

        found_split_guides = []
        for offset in offsets:
            config = self._create_split_config(offset)
            split_guides = guide_splitter.split(config)
            if self._has_all_unique_overlaps(split_guides):
                found_split_guides = split_guides
                break
            else:
                continue

        return found_split_guides

    def _create_split_config(self, offset):
        centers = [guide.center for guide in self.guides]
        return SplitConfiguration([sum(x) for x in zip(offset, centers)])

    def _has_all_unique_overlaps(self, split_guides):
        overlaps = [split_guide.overlap for split_guide in split_guides]
        overlaps_set = set()
        for overlap in overlaps:
            overlaps_set.add(str(overlap))
            overlaps_set.add(str(overlap.reverse_complement()))

        if len(overlaps_set) == 2 * len(split_guides):
            return True
        return False
