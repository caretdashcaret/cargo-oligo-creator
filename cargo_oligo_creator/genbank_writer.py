from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

class GenbankWriter:

    def __init__(self, path, displays):
        self.displays = displays
        self.path = path
        self.index = 0
        self.success = False

    def write(self):
        for display in self.displays:
            self._write_display(display)
        self.success = True

    def _write_display(self, display):
        self._write_oligo(display, display.forward_oligo)
        self._write_oligo(display, display.reverse_oligo)

    def _write_oligo(self, display, oligo):
        sequence = Seq(oligo.raw_sequence(), IUPAC.unambiguous_dna)

        record = SeqRecord(sequence,
                           id=str(self.index),
                           name=self._get_record_name(display, oligo),
                           description=self._get_description(display, oligo))

        filename_path = self._get_full_filepath(self._get_filename(display, oligo))

        output_file = open(filename_path,"w")
        SeqIO.write(record, output_file, "genbank")
        self.index += 1

    def _get_filename(self, display, oligo):
        return "{}_oligo_for_{}_to_{}.gb".format(oligo.direction(),
                                                 display.pre_construct_guide_index,
                                                 display.post_construct_guide_index)

    def _get_record_name(self, display, oligo):
        return "{}_to_{}:{}".format(display.pre_construct_guide_index, display.post_construct_guide_index, oligo.direction())

    def _get_description(self, display, oligo):
        return "{}:{}".format(display.name(), oligo.display_formatted())

    def _get_full_filepath(self, filename):
        return os.path.join(os.path.abspath(self.path), filename)
