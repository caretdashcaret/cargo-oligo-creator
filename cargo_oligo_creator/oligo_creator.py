from Bio.Seq import Seq
from .oligo import Oligo
from .oligos_display import OligosDisplay

class OligoCreator:
    # The BsaI constructs sit inbetween the 2nd part of the n-th gRNA and the 1st part of the n+1th gRNA
    #
    # 5' --- NNNNNGAGACC|NN|GGTCTCN     --- 3'
    # 3' ---     NCTCTGG|NN|CCAGAGNNNNN --- 5'

    # Here we define our construct to be:

    # 5' --- NNNNTGAGACC|GA|GGTCTCA     --- 3'
    # 3' ---     ACTCTGG|CT|CCAGAGTNNNN --- 5'

    BSAI_CONSTRUCT_FORWARD = Seq("T") + Seq("GAGACC") + Seq("GA") + Seq("GGTCTC") + Seq("A")
    BSAI_CONSTRUCT_REVERSE = BSAI_CONSTRUCT_FORWARD.reverse_complement()

    FORWARD_START = Seq("CACC")
    REVERSE_START = Seq("AAAC")

    def __init__(self, split_guides):
        self.split_guides = split_guides

    def create_oligos(self):
        if len(self.split_guides) == 0:
            return []
        else:
            return self._join_to_form_oligos(self.split_guides)

    def _create_forward_and_reverse_oligos(self, split_guide_a, split_guide_b):
        forward_oligo = self._create_forward_oligo(split_guide_a, split_guide_b)
        reverse_oligo = self._create_reverse_oligo(split_guide_a, split_guide_b)
        return [forward_oligo, reverse_oligo]

    def _create_forward_oligo(self, split_guide_a, split_guide_b):
        return Oligo(
            is_forward = True,
            overlap = split_guide_b.overlap,
            start = self.FORWARD_START,
            pre_construct_piece = split_guide_a.first_part,
            construct = self.BSAI_CONSTRUCT_FORWARD,
            post_construct_piece = split_guide_b.second_part)

    def _create_reverse_oligo(self, split_guide_a, split_guide_b):
        return Oligo(
            is_forward = False,
            overlap = split_guide_b.overlap.reverse_complement(),
            start = self.REVERSE_START,
            pre_construct_piece = split_guide_b.second_part.reverse_complement(),
            construct = self.BSAI_CONSTRUCT_REVERSE,
            post_construct_piece = split_guide_a.first_part.reverse_complement())

    def _join_to_form_oligos(self, split_guides):
        oligos = []
        index = 0
        while index < len(split_guides):
            post_construct_guide_index = self._get_index_of_post_construct_guide(index)
            forward_oligo, reverse_oligo = self._create_forward_and_reverse_oligos(self.split_guides[index],
                                                                                   self.split_guides[post_construct_guide_index])
            oligos.append(OligosDisplay(
                forward_oligo = forward_oligo,
                reverse_oligo = reverse_oligo,
                pre_construct_guide_index = index,
                post_construct_guide_index = post_construct_guide_index))

            index += 1
        return oligos

    def _is_first_guide(self, guide_index):
        return guide_index == 0

    def _get_index_of_post_construct_guide(self, index):
        # if it's the 1st guide, the part after BSAI should come from the last guide
        if index == 0:
            return len(self.split_guides) - 1
        else:
            return index - 1
