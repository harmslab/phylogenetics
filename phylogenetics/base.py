#Base classes for Homolog objects
#
#

import numpy as np
import json
import pickle


class Homolog(object):

    def __init__(self, unique_id, **kwargs):
        """ Create a single Homolog object.

            Args:
            ----
            unique_id: str
                ID number for homolog, unique from any set.
        """
        # Must set a unique ID and sequence
        self.id = unique_id

        # Set user specified attributes
        for key, value in kwargs.items():
            # Protect from overwriting Python native attributes
            if key == "def":
                setattr(self, "defline", value)
            elif key == "id":
                setattr(self, "gid", value)
            else:
                setattr(self, key, value)

    def add_attributes(self, **kwargs):
        """ Add attributes to homolog object. """
        for key,value in kwargs.items():
            setattr(self, key, value)

    # -----------------------------------
    # Output formats
    # -----------------------------------

    def fasta(self, tags=None, aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        f = ">"
        if tags is not None:
            # Remove sequence if its an elements
            if "sequence" in tags:
                tags.remove("sequence")
            f += "|".join([str(getattr(self, t)) for t in tags])
        else:
            f += self.id

        # New line with full sequence string (aligned if told to)
        if aligned:
            f += "\n" + self.latest_align +"\n"
        else:
            f += "\n" + self.sequence +"\n"

        return f

    def phylip(self, tags=None, **kwargs):
        """ Return a PhyML formatted string to write to file.

            Only allowed when Homolog has alignment attribute
        """
        # Check that an alignment attribute exists in homolog.
        if hasattr(self, "latest_align") is False:
            raise Exception("Homolog must have an alignment attribute to write phylip.")

        # Use specified alignment key
        if tags is not None:
            alignment = getattr(self, tags)
        else:
            alignment = self.latest_align

        f = "%s\n%s\n" % (self.id, alignment)
        return f

    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.dumps(self.__dict__)

    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self)

    def csv(self, tags=None, header=True, delimiter=",", **kwargs):
        """ write csv string. """
        # Get all attributes if tags are not specified
        if tags is None:
            tags = list(vars(self).keys())

        # Add header is specified
        if header:
            f = delimiter.join(tags)
        else:
            f = ""

        # Try to print tag. If it doesnt exist, then leave blank
        vals = list()
        for t in tags:
            try:
                vals.append(str(getattr(self,t)))
            except:
                vals.append("")
        f += delimiter.join(vals) + "\n"

        return f

    def write(self, filename, format="fasta", tags=None, aligned=False):
        """ Write to file with given format.

            Default is fasta.
        """
        # Different writing formats that are possible
        write_format = {"fasta":"w",
                        "pickle": "wb",
                        "json":"w",
                        "phylip":"w"}

        format_func = getattr(self, format)
        f = open(filename, write_format[format])
        f.write(format_func(tags=tags, aligned=aligned))
        f.close()


class HomologSet(object):

    def __init__(self, homolog_set=[]):
        """ Construct a set of homolog objects. """
        self._homologs = homolog_set

    @property
    def homologs(self):
        """ Get homolog set. """
        return self._homologs

    def get_map(self, attr1, attr2=None):
        """ Return mapping between two attributes in homolog set, OR
            if no second attribute is give, map attr1 to whole homolog.

            attr2 can be a list of other attributes -- which will return
                a list of those attributes mapped to attr1.
        """
        m = dict()

        # If no second attribute is give, mapping is between first attribute
        # and the homolog object
        if attr2 is None:
            for h in self._homologs:
                m[getattr(h, attr1)] = h

        # else, mapping from one attribute to another
        else:

            # If attr2 is a list, return a list in mapping
            if isinstance(attr2,list):
                for h in self._homologs:
                    m[getattr(h, attr1)] = [getattr(h, a) for a in attr2]

            # Else just return a single attr.
            else:
                for h in self._homologs:
                    m[getattr(h, attr1)] = getattr(h, attr2)

        return m

    def add_homologs(self, homologs):
        """ Append a list of homolog objects to the set."""

        # If a single homolog is given, format it into a list
        # for loop below.
        if isinstance(homologs,list) == False:
            homologs = [homologs]

        # Append each homolog to list if it's a homolog instance.
        for i in range(len(homologs)):
            if isinstance(homologs[i], Homolog):
                self._homologs.append(homologs[i])
            else:
                raise Exception("homolog must be an instance of Homolog class.")

    def rm_homologs(self, ids):
        """ Remove a list of homologs from set of homologs."""
        # If a single id is given, format it into a list
        # for loop below.
        if isinstance(ids,list) == False:
            ids = [ids]

        # Get id map
        homolog_dict = self.get_map("id")

        # Remove these homologs from system
        for id in ids:
            del homolog_dict[id]

        self._homologs = list(homolog_dict.values())


    def subset_homologs(self, ids):
        """ Returns a new Homologs object from subset of this homolog object. """
        if isinstance(ids,list) == False:
            ids = [ids]

        # Get id map
        homolog_dict = self.get_map("id")

        # Remove these homologs from system
        subset = list()
        for id in ids:
            subset.append(homolog_dict[id])

        return HomologSet(homolog_set=subset)

    def find_homolog(self, id):
        """ Returns the metadata dictionary of homolog with id. """
        mapping = self.get_map("id")
        return mapping[id].__dict__

    def print_homolog(self, id, keys=None):
        """ Print the metadata dictionary of homolog with id. """

        # Homolog dictionary
        homolog_dict = self.find_homolog(id)
        if keys is None:
            keys = list(homolog_dict.keys())

        # Print in table
        print(id + ":\n" + "-----------\n")
        for key in keys:
            print(key+" : " + str(homolog_dict[key]) + "\n")

    # -----------------------------------
    # Output formats
    # -----------------------------------

    def fasta(self, tags=None, aligned=False):
        """ Return string in fasta format for the set."""
        f = ""
        for h in self._homologs:
            f += h.fasta(tags, aligned=aligned)
        return f

    def phylip(self, tags=None, **kwargs):
        """ Return string of sequences in phylip format. """

        # Get the latest align if other alignment isn't specified
        if tags is None:
            tags = "latest_align"

        f = ""
        for h in self.homologs:
            f += h.phylip(tags)

        n_homologs = len(self.homologs)
        n_col = len(getattr(self.homologs[0], tags))

        out = "%i  %i\n\n%s\n" % (n_homologs,n_col,f)
        return out

    def json(self, **kwargs):
        """ Return json string of homolog set."""
        obj = list()
        for h in self._homologs:
            obj.append(h.__dict__)
        return json.dumps(obj)

    def pickle(self, **kwargs):
        """ Return pickle string of homolog set. """
        return pickle.dumps(self)

    def csv(self, tags=None, delimiter=",", **kwargs):
        """ Return csv string. """
        # If tags is not specified, get all tags.
        if tags is None:
            tags = list(self.homologs[0].__dict__.keys())

        f = delimiter.join(tags)
        f += "\n"
        for h in self._homologs:
            f += h.csv(tags=tags, header=False, delimiter=delimiter)
        return f

    def write(self, filename, format="fasta", tags=None, aligned=False):
        """ Write to file with given format.

            Default is fasta.
        """
        # Different writing formats that are possible
        write_format = {"fasta":"w",
                        "pickle": "wb",
                        "json":"w",
                        "phylip":"w",
                        "csv": "w"}

        format_func = getattr(self, format)
        f = open(filename, write_format[format])
        f.write(format_func(tags=tags, aligned=aligned))
        f.close()

# -----------------------------------------------------
# Functions to manage and maintain homologs.
# -----------------------------------------------------

def rank_homologs(homolog_set, positive=(), negative=("putative","hypothetical","unnamed",
                    "possible", "predicted","unknown","uncharacterized",
                    "mutant","isoform"), rank_offset=0):

    """ Rank homologs based on dubious descriptions in their defline. """

    for h in homolog_set.homologs:
        defline = h.defline

        # Does one of the dubious entries occur on this line?
        rank = rank_offset

        for p in positive:
            # If positive strings are in defline, subtract from rank
            if p in defline:
                rank -= 1

        for n in negative:
            # If negative strings are in defline, add to rank
            if n in defline:
                rank += 1

        h.add_attributes(rank=rank)
