class Oligo(object):

    def __init__(self, is_forward, overlap, guide_end, guide_start, start, pre_construct_piece, construct, post_construct_piece):
        self.is_forward = is_forward
        self.overlap = overlap
        self.guide_end = guide_end
        self.guide_start = guide_start
        self.start = start
        self.pre_construct_piece = pre_construct_piece
        self.construct = construct
        self.post_construct_piece = post_construct_piece
        self.formatted_sequence = self.start + "|" + self.pre_construct_piece + "|" + self.construct + "|" + self.post_construct_piece


    def __repr__(self):
        direction = "Forward"
        if not self.is_forward:
            direction = "Reverse"

        return "{} oligo for [end of guide {}] to [start of guide {}] with overlap {}: {}".format(
            direction,
            self.guide_end,
            self.guide_start,
            self.overlap,
            self.formatted_sequence)
