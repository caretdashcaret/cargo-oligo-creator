from Bio.Seq import Seq
from .oligo import Oligo
from .split_guides_finder import SplitGuidesFinder

class OligoCreator:
    # The BsaI constructs sit inbetween the 2nd part of the n-th gRNA and the 1st part of the n+1th gRNA
    #
    # 5' --- NNNNNGAGACC|NN|GGTCTCN     --- 3'
    # 3' ---     NCTCTGG|NN|CCAGAGNNNNN --- 5'

    # Here we define our construct to be:

    # 5' --- NNNNTGAGACC|GA|GGTCTCA     --- 3'
    # 3' ---     ACTCTGG|CT|CCAGAGTNNNN --- 5'

    BSAI_CONSTRUCT_FORWARD = Seq("T") + Seq("GAGACC") + Seq("GA") + Seq("GGTCTC") + Seq("A")
    BSAI_CONSTRUCT_REVERSE = Seq("A") + Seq("CTCTGG") + Seq("CT") + Seq("CCAGAG") + Seq("T")

    FORWARD_START = Seq("CACC")
    REVERSE_START = Seq("AAAC")

    def __init__(self, guides):
        self.guides = guides

    def create_oligos(self):
        found_split_guides = SplitGuidesFinder(self.guides).find_split_guides()

        if len(found_split_guides) == 0:
            return []
        else:
            return self._join_to_form_oligos(found_split_guides)

    def _join_to_form_oligos(self, split_guides):
        oligos = []
        index = 0
        while index < len(split_guides):
            current = split_guides[index]
            if not self._is_first_guide(index):
                previous = split_guides[index-1]
                forward_oligo = Oligo(
                    is_forward = True,
                    overlap = current.overlap,
                    guide_end = index - 1,
                    guide_start = index,
                    start = self.FORWARD_START,
                    pre_construct_piece = previous.second_part,
                    construct = self.BSAI_CONSTRUCT_FORWARD,
                    post_construct_piece = current.first_part,
                )
                reverse_oligo = Oligo(
                    is_forward = False,
                    overlap = previous.second_part[-4:],
                    guide_end = index - 1,
                    guide_start = index,
                    start = self.REVERSE_START,
                    pre_construct_piece = current.first_part.reverse_complement(),
                    construct = self.BSAI_CONSTRUCT_REVERSE,
                    post_construct_piece = previous.second_part.reverse_complement()
                )
                oligos.append(forward_oligo)
                oligos.append(reverse_oligo)
            index += 1
        return oligos


    def _is_first_guide(self, guide_index):
        return guide_index == 0
