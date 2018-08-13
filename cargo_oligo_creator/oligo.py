class Oligo(object):

    def __init__(self, is_forward, overlap, start, pre_construct_piece, construct, post_construct_piece):
        self.is_forward = is_forward
        self.overlap = overlap
        self.start = start
        self.pre_construct_piece = pre_construct_piece
        self.construct = construct
        self.post_construct_piece = post_construct_piece

    def formatted_sequence(self):
        return self.start + "|" + self.pre_construct_piece + "|" + self.construct + "|" + self.post_construct_piece

    def __repr__(self):
        direction = "Forward"
        if not self.is_forward:
            direction = "Reverse"

        return "{} oligo with overlap {}: {}".format(
            direction,
            self.overlap,
            self.formatted_sequence())
