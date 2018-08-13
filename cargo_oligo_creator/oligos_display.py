class OligosDisplay:

    # This should hold the same information as oligo except for additional information about the original guides
    def __init__(self, forward_oligo, reverse_oligo, pre_construct_guide_index, post_construct_guide_index):
        self.forward_oligo = forward_oligo
        self.reverse_oligo = reverse_oligo
        self.pre_construct_guide_index = pre_construct_guide_index
        self.post_construct_guide_index = post_construct_guide_index

    def name(self):
        return "Oligos for [start piece of guide {}] to [end piece of guide {}]".format(
            self.pre_construct_guide_index,
            self.post_construct_guide_index
        )
