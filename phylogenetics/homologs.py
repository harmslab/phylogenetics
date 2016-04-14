# API for working with Homolog sets in a phylogenetics project

import json
import pickle
from phylogenetics.names import switch

from phylogenetics.dataio.write import Write, WriteSet

# ---------------------------------------------------
# Things that you often do with HomologSets
# ---------------------------------------------------

def unique_id(start, end=None):
    """ Returns a set of unique names """
    if end is None:
        end = int(start)
        start = 0
    return ["XX%08d" % i for i in range(start, end+1)]

def load_homologset(filename):
    """ Load a homologset from pickle file.

        Note: only need to give the filename.
    """
    f = open(filename, "rb")
    homologset = pickle.load(f)
    f.close()
    return homologset

def concat_homolog_sets(hs1, hs2, renumber=False):
    """ Concatenate two homolog set. If told, will renumber the `id` attributes
        in the merged set.
    """
    total_set = hs1.homologs + hs2.homologs
    total_hs = HomologSet(total_set)
    if renumber:
        total_hs.renumber_homologs()
    return total_hs

def rm_repeats_homologs(homolog_set, attribute="accession", renumber=False):
    """ Remove and repetitive homologs in a set with
        respect to a given attribute.

        Arguments:
        ---------
        homolog_set: HomologSet object
            Set to search through.
        attribute: str (default="accession")
            Attribute for using to search repeats in set

    """
    # Get homologs from set
    homs = homolog_set.homologs
    # Get the attibutes in the order of the homologs
    attributes = [getattr(h,attribute) for h in homs]

    seen = set() # a set for keeping already seen attributes
    unique_homologs = [] # The list for storing homolog subset

    # Iterate through list and find unique homologs
    for i in range(len(homs)):
        if attributes[i] not in seen:
            unique_homologs.append(homs[i])
            seen.add(attributes[i])

    # Build a new set of homologs
    hs = HomologSet(unique_homologs)

    # renumber ID's if told to.
    if renumber:
        hs.renumber_homologs()
    return hs


def rank_homologs(homolog_set, accession=(), positive=(), negative=("putative","hypothetical","unnamed",
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

        # If accession is given, it should supercede all other ranks.
        try:
            access = h.accession
            for a in accession:
                if a in access:
                    rank -= 100
        except:
            pass

        h.add_attributes(rank=rank)

# ---------------------------------------------------
# Main Homolog objects for package
# ---------------------------------------------------

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
        self.latest_align = ""

        # Attach Write-ing object
        self.Write = Write(self)

        self.__tags = {}
        # Set user specified attributes
        for key, value in kwargs.items():
            # Protect from overwriting Python native attributes
            if key == "def":
                self.__tags["defline"] = value
            elif key == "id":
                self.__tags["defline"] = value
            else:
                self.__tags[key] = value


    @property
    def __tags__(self):
        """ """
        return self.__tags

    @property
    def seqlen(self):
        """ Get the length of sequence """
        return len(self.sequence)

    @property
    def alignedlen(self):
        """ Get the length of the aligned sequence. """
        return len(self.latest_align) - self.latest_align.count('-')

    def add_attributes(self, **kwargs):
        """ Add attributes to homolog object. """
        for key,value in kwargs.items():
            setattr(self, key, value)


class HomologSet(object):

    def __init__(self, homolog_set=[]):
        """ Construct a set of homolog objects. """
        self._homologs = homolog_set

        # Attach Write-ing object
        self.Write = WriteSet(self)

    @property
    def homologs(self):
        """ Get homolog set. """
        return self._homologs

    @property
    def id_list(self):
        """ Return ID list. """
        return [h.id for h in self._homologs]

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

    def subset(self, ids):
        """ Get a subset of the homologs from an ID list. """
        mapping = self.get_map("id")

        homologs = []
        for id in ids:
            homologs.append(mapping[id])

        return HomologSet(homologs)

    def renumber_homologs(self):
        """ Renumber the ID numbers for homologs in the set, starting at
            0 to len(homologs).
        """
        n = len(self.homologs)
        for i in range(n):
            unique_id = "XX%08d" % i
            self._homologs[i].id = unique_id

    def add_homologs(self, homologs):
        """ Append a list of homolog objects to the set.

            NOTE: does not renumber the homolog set. this must be called
            manually.
        """

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
        """ Remove a list of homologs from set of homologs.

            NOTE: this does not renumber the homolog set. must be done manually.

            Also, note that this changes the homolog set in place!
        """
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
        return mapping[id].__tags__

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
