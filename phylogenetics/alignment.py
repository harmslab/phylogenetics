## Tools for managing alignments

from phylogenetics.dataio.alignmentio import Write, Read

class Alignment(object):

    def __init__(self, HomologSet):
        """ Object for maintaining alignment data for a HomologSet Object.
        """
        self._HomologSet = HomologSet
        self.Write = Write(self)
        self.Read = Read(self)

    @property
    def latest(self):
        """ Return the latest alignment. """
        self.Write.fasta()

    @property
    def n_taxa(self):
        """ Return the number of sequences in alignment """
        return len(self.sequences)

    @property
    def length(self):
        """ Return the length of the alignment. """
        return len(self.sequences[0])

    @property
    def without_gaps(self):
        """ """
        pass

    def update(self, alignment="latest_align"):
        """ Update an alignment with new set of sequences. """
        pass
