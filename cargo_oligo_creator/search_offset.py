class SearchOffset:

    # Search at most SEARCH_LIMIT positions from the inital starting position for unique overlaps (which become overhangs)
    # This can be modified but the running time is O( (2+1+SEARCH_LIMIT)^N ) where N is is the number of guides
    SEARCH_LIMIT = 3

    def __init__(self, max_number_of_guides):
        self.max_number_of_guides = max_number_of_guides
        self._ordered_offsets = self._generate_ordered_offsets_for_single_guide_search()

    def get_offsets(self):
        return self.generate_offsets_for_search(0)

    def _generate_ordered_offsets_for_single_guide_search(self):
        # For a SEARCH_LIMIT of 3, this generates [0, -1, 1, -2, 2, -3, 3]
        offsets = list(range(-1 * self.SEARCH_LIMIT, self.SEARCH_LIMIT + 1))
        return sorted(offsets, key=abs)

    def generate_offsets_for_search(self, guide_index):
        if self._is_last_guide(guide_index):
            return [[x] for x in self._ordered_offsets]
        else:
            offsets = []
            current_offsets = [[x] for x in self._ordered_offsets]
            rest_of_offsets = self.generate_offsets_for_search(guide_index+1)
            for each_offset in current_offsets:
                for rest_of_offset in rest_of_offsets:
                    offsets.append(each_offset + rest_of_offset)
            return offsets

    def _is_last_guide(self, guide_index):
        return guide_index == self.max_number_of_guides - 1
