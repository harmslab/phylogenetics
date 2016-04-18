# Module for input and output of alignment.

from .base import read_from_file, write_to_file
from . import fasta

class Write(object):

    def __init__(self, Alignment):

        self._Alignment = Alignment

    def fasta(self, name="latest_align"):
        """ Write alignment to fasta. """
        pass

    def phylip(self):
        pass

    def csv(self):
        pass
