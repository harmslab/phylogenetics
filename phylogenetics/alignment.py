## Tools for managing alignments

from phylogenetics.dataio.alignmentio import Write, Read

class Alignment(object):

    def __init__(self, HomologSet):
        """Object for maintaining alignment data for a HomologSet Object.
        """
        self._HomologSet = HomologSet
        self.Write = Write(self)
        self.Read = Read(self)
        self._alignments = {"latest":{}}

    @property
    def latest(self):
        """ Return the latest alignment. """
        return self._alignments["latest"]

    @property
    def n_taxa(self):
        """ Return the number of sequences in alignment """
        return len(self.latest)

    @property
    def length(self):
        """ Return the length of the alignment. """
        return len(list(self.latest.values())[0])

    def list_alignments(self):
        """List the key for all alignments stored in this object."""
        return list(self._alignments.keys())

    def previous(self, key):
        """Get previous alignments
        """
        return self._alignments[key]

    def _move_old_alignment(self):
        """Moves the current alignment to new attribute name to make room for
        a alignment.
        """
        if self.latest != {}:
            # Save old alignments.
            exists = True
            counter = 0
            while exists == True:
                if "align" + str(counter) in self._alignments:
                    counter += 1
                else:
                    exists = False

            # Move old alignment to new name.
            self._alignments["align" + str(counter)] = self.latest

    def update(self, alignment_file):
        """ Update an alignment with new set of sequences.

        Reads in an alignment file and moves old alignments to new attribute.
        """
        # Write out alignment file
        self.Read.fasta(fname=alignment_file)
