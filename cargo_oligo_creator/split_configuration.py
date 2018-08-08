class SplitConfiguration:

    def __init__(self, split_indexes):
        self.split_indexes = split_indexes

    def get_split_index_for_guide(self, guide_index):
        return self.split_indexes[guide_index]
