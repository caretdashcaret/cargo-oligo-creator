class Oligo(object):

    def __init__(self, is_forward, overlap, start, pre_construct_piece, construct, post_construct_piece):
        self.is_forward = is_forward
        self.overlap = overlap
        self.start = start
        self.pre_construct_piece = pre_construct_piece
        self.construct = construct
        self.post_construct_piece = post_construct_piece

    def _sequence_components(self):
        return [str(x) for x in [self.start, self.pre_construct_piece, self.construct, self.post_construct_piece]]

    def formatted_sequence(self):
        return "|".join(self._sequence_components())

    def raw_sequence(self):
        return "".join(self._sequence_components())

    def display_formatted(self):
        return "{} oligo with overlap {}: {}".format(
            self._direction(),
            self.overlap,
            self.formatted_sequence())

    def display_raw(self):
        return "{} oligo with overlap {}: {}".format(
            self._direction(),
            self.overlap,
            self.raw_sequence())

    def _direction(self):
        if self.is_forward:
            return "Forward"
        else:
            return "Reverse"

    def __repr__(self):
        return self.display_formatted()
