## Tools for managing alignments

class Alignment(object):

    def __init__(self, HomologSet):
        """ Object for maintaining alignment data for a HomologSet Object.
        """
        self._HomologSet = HomologSet

    @property
    def latest(self, fname=None):
        """ Return the latest alignment. """
        self._homologset.write.fasta(fname=fname, aligned=True)

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
