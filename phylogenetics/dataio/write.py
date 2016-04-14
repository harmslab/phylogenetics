# Module for writing HomologSets out to various formats

import json
import pickle

from functools import wraps

def _write(function):
    """ Decorator for determining whether to write to file or print to terminal."""

    #@wraps # functools for ridding ourselves of closure for pickling
    def wrapper(self, fname=None, *args, **kwargs):
        """ """
        # Create a string of whatever datatype
        string = function(self, *args, **kwargs)

        # If a filename is not given, return string.
        if fname is None:
            return string

        # Write to a file
        else:
            try:
                # Try to write a straight string.
                with open(fname, "w") as f:
                    f.write(string)
            except:
                # IF writing string failed, try writing bytes.
                with open(fname, "wb") as f:
                    f.write(string)

            return string

    return wrapper


class Write(object):

    def __init__(self, Homolog):
        """ """
        self._homolog = Homolog

    @_write
    def fasta(self, tags=None, aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        f = ">"
        if tags is not None:
            # Remove sequence if its an elements
            if "sequence" in tags:
                tags.remove("sequence")
            f += "|".join([str(getattr(self._homolog, t)) for t in tags])
        else:
            f += self.id

        # New line with full sequence string (aligned if told to)
        if aligned:
            f += "\n" + self._homolog.latest_align +"\n"
        else:
            f += "\n" + self._homolog.sequence +"\n"

        return f

    @_write
    def phylip(self, tags=None, **kwargs):
        """ Return a PhyML formatted string to write to file.

            Only allowed when Homolog has alignment attribute
        """
        # Check that an alignment attribute exists in homolog.
        if hasattr(self._homolog, "latest_align") is False:
            raise Exception("Homolog must have an alignment attribute to write phylip.")

        # Use specified alignment key
        if tags is not None:
            alignment = getattr(self._homolog, tags)
        else:
            alignment = self._homolog.latest_align

        f = "%s\n%s\n" % (self._homolog, alignment)
        return f


    @_write
    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.dumps(self._homolog.__tags__)


    @_write
    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self._homolog)


    @_write
    def csv(self, tags=None, header=True, delimiter=",", **kwargs):
        """ write csv string. """
        # Get all attributes if tags are not specified
        if tags is None:
            tags = list(vars(self._homolog).keys())

        # Add header is specified
        if header:
            f = delimiter.join(tags)
        else:
            f = ""

        # Try to print tag. If it doesnt exist, then leave blank
        vals = list()
        for t in tags:
            try:
                vals.append(str(getattr(self._homolog,t)))
            except:
                vals.append("")
        f += delimiter.join(vals) + "\n"

        return f



class WriteSet(object):

    def __init__(self, HomologSet):
        """ Object for writing out metadata held in a HomologSet. """
        self._homologset = HomologSet

    @_write
    def fasta(self, tags=None, aligned=False):
        """ Return string in fasta format for the set."""
        f = ""
        for h in self._homologset._homologs:
            f += h.Write.fasta(tags, aligned=aligned)
        return f

    @_write
    def phylip(self, tags=None, **kwargs):
        """ Return string of sequences in phylip format. """

        # Get the latest align if other alignment isn't specified
        if tags is None:
            tags = "latest_align"

        f = ""
        for h in self._homologset.homologs:
            f += h.Write.phylip(tags)

        n_homologs = len(self._homologset.homologs)
        n_col = len(getattr(self._homologset.homologs[0], tags))

        out = "%i  %i\n\n%s\n" % (n_homologs,n_col,f)
        return out

    @_write
    def json(self, **kwargs):
        """ Return json string of homolog set."""
        obj = list()
        for h in self._homologset._homologs:
            obj.append(h.__tags__)
        return json.dumps(obj)

    @_write
    def pickle(self, **kwargs):
        """ Return pickle string of homolog set. """
        return pickle.dumps(self._homologset)

    @_write
    def csv(self, tags=None, delimiter=",", **kwargs):
        """ Return csv string. """
        # If tags is not specified, get all tags.
        if tags is None:
            tags = list(self._homologset.homologs[0].__tags__.keys())

        f = delimiter.join(tags)
        f += "\n"
        for h in self._homologset._homologs:
            f += h.Write.csv(tags=tags, header=False, delimiter=delimiter)
        return f

    @_write
    def newick(self, tags, **kwargs):
        """ Write a tree to file. """
        old_name = "id"
        new_names = tags
        tree = switch(self._homologset, "id", new_names, format="newick")
        return tree
