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
        return self.Write.fasta()

    @property
    def n_taxa(self):
        """ Return the number of sequences in alignment """
        return len(self._HomologSet.homologs)

    @property
    def length(self):
        """ Return the length of the alignment. """
        return len(list(self._HomologSet.homologs.values())[0].latest_align)

    def update(self, alignment_file):
        """ Update an alignment with new set of sequences. """
        # Write out alignment file
        self.Read.fasta(fname=alignment_file)
