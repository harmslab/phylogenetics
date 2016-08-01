# API for working with Homolog sets in a phylogenetics project
from __future__ import absolute_import

import os, string

# Reading and Writing modules to HomologSet
import phylogenetics.dataio.homologio as homologio
import phylogenetics.dataio.homologsetio as homologsetio

# Subobjects bound to HomologSet
from phylogenetics.alignment import Alignment
from phylogenetics.tree import Tree
from phylogenetics.ancestors import Ancestor, AncestorSet

# Import external tools
from phylogenetics.dataio.formats import entrez_xml
from phylogenetics.exttools import (cdhit,
                                msaprobs,
                                phyml,
                                paml,
                                entrez,
                                taxonomy)

import phylogenetics.utils as utils

# ---------------------------------------------------
# Things that you often do with HomologSets
# ---------------------------------------------------

def unique_id(start, end=None):
    """ Returns a set of unique names """
    if end is None:
        end = int(start)
        start = 0
    return ["XX%08d" % i for i in range(start, end+1)]

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

        Parameters
        ----------
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


def rank_homologs(homolog_set,
        accession=(),
        positive=(),
        negative=(
            "putative","hypothetical","unnamed","possible", "predicted",
            "unknown","uncharacterized","mutant","isoform"
        ),
        rank_offset=0
    ):
    """Rank homologs based on dubious descriptions in their defline.
    """
    for id, h in homolog_set.homologs.items():
        h.addattr("rank", rank_offset)
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
            if hasattr(h, "accver"):
                access = h.accver
            else:
                access = h.accession
            for a in accession:
                if a in access:
                    rank -= 100
        except:
            pass

        h.rank = rank

# ---------------------------------------------------
# Main Homolog objects for package
# ---------------------------------------------------

class Homolog(object):

    def __init__(self, unique_id, **kwargs):
        """ Inialize a Homolog Object.

            The `Homolog` object provides a data-structure for sequence metadata
            in a phylogenetics project.


            Arguments:
            ---------
            unique_id: str
                ID number for homolog, unique from any set.

            kwargs become tags/
        """
        # Initialize all attributes
        self._attrs = {}

        # Must set a unique ID and sequence
        self.addattr("id", unique_id)

        # Set user specified attributes
        for key, value in kwargs.items():
            # Protect from overwriting Python native attributes
            if key == "def":
                self.addattr("defline",value)
            elif key == "id":
                self.attattr("gid", value)
            else:
                self.addattr(key, value)

        # Attach Write-ing object
        self.Write = homologio.Write(self)
        self.Read = homologio.Read(self)

    @classmethod
    def download(self, accession):
        """ """

    @property
    def attrs(self):
        """ """
        return self._attrs

    @property
    def seqlen(self):
        """ Get the length of sequence """
        return len(self.sequence)

    def get(self, *attrs):
        """ Get specific attributes from Homolog.

        Arguments
        ---------
        *attrs : list
            list of attributes to retrieve from Homolog object.

        If no attrs are given, return all attrs.

        Returns
        -------
        metadata : dict
            returns attributes that were given.
        """
        metadata = {}
        for a in attrs:
            metadata[a] = self._attrs[a]

        if len(metadata) == 0:
            metadata = self.attrs

        return metadata

    def addattr(self, key, value):
        """ Add attributes to homolog object. """
        # Set in attrs dict
        self._attrs[key] = value

        # Set the object attribute as well
        setattr(self, key, value)

    def rmattr(self, attr):
        """ remove an attribute from Homologs. """
        # Remove from _attrs dict
        del self._attrs[attr]

        # Remove from object
        delattr(self, attr)


class HomologSet(object):

    def __init__(self, homologs=[]):
        """ Initialize a HomologSet object.

            A HomologSet object provides a data-structure for metadata from a set
            of sequences for a phylogenetics project.

        """
        self._homologs = {}
        self.add(homologs)

        # Attach Write-ing object
        self.Write = homologsetio.Write(self)
        self.Read = homologsetio.Read(self)

        setattr(self, "download", self._download)

    @classmethod
    def download(cls, accessions, email):
        """ Download a HomologSet from Entrez
        """
        # Check that ids is a list
        if type(accessions) != list:
            raise Exception("""`accessions` must be a list.""")

        # Initialize the class
        hs = cls()
        # Use download method inside class
        hs._download(accessions, email)
        # overwrite this classmethod with non-classmethod
        hs.download = hs._download
        return hs

    def _download(self, accessions, email):
        """ Download HomologSet
        """
        # Check that ids is a list
        if type(accessions) != list:
            raise Exception("""`accessions` must be a list.""")

        # Download the full metadata for
        data = entrez.download(accessions, email)

        self.Read.entrez_xml(data)

        # Quality control for data without sequences in them.
        bad_data = [id for id, Homolog in self.homologs.items() if hasattr(Homolog, "sequence") is False]

        # Remove bad data from HomologSet
        self.rm(bad_data)

    def download_taxonomy(self, greedy=False):
        """Retrieve Taxonomic information about HomologSet.

        Note: This removes sequences that do not have taxonomic data.
        """
        first = True
        bad_data = []
        # Download taxonomy data from BLAST
        for h in self.homologs:
            # If sequence doesn't have an organism name, remove
            homolog = self.homologs[h]
            try:
                # Estimate how long this will take
                if first:
                    time = utils.timeit(taxonomy.query, homolog.taxid)
                    print("This will take roughly %d seconds" % round(time*self.length))

                output = taxonomy.query(homolog.taxid)
                homolog.addattr("taxonomy", output)
                first = False
            except AttributeError:
                bad_data.append(homolog.id)

        print("Done.")

        if greedy is True:
            # Remove bad data from HomologSet
            self.rm(bad_data)

    @classmethod
    def load(cls, fname):
        """ Read a HomologSet. """
        with open(fname, "rb") as f:
            homologset = pickle.load(f)
        return homologset

    def retrieve(self, ids, email):
        """Retrieve data from Entrez for Homologs already
        in the HomologSet object.
        """
        # Check that ids is a list
        if type(ids) != list:
            raise Exception("""`ids` must be a list.""")

        # Mapping to accession
        mapping = self.map("id", "accver")

        # List accessions.
        accessions = [mapping[id] for id in ids]
        self._download(accessions, email)

    @property
    def length(self):
        """Return the number of homologs in the HomologSet."""
        return len(self._homologs)

    @property
    def homologs(self):
        """ Get homolog set. """
        return self._homologs

    @property
    def list_ids(self):
        """ Return ID list. """
        return [hom.id for h, hom in self._homologs.items()]

    @property
    def max_id(self):
        """ Return the max id. """
        if len(self.list_ids) == 0:
            return 0
        else:
            return max([int(id[2:]) for id in self.list_ids])

    def get(self, *attrs):
        """Get a dictionary with metadata for all homologs objects in the
        HomologSet. The metadat for each Homolog only includes the attributes
        named as arguments.

        Arguments
        ---------
        *attrs : strings
            Arguments to include in dictionary

        If no attributes are given, return all.

        Returns
        -------
        metadata : dict
            Dictionary will metadata for each Homolog in HomologSet; only includes
            keys named in arguments.
        """
        metadata = {}
        for h in self.list_ids:
            # Get homolog object
            Homolog = getattr(self, h)
            metadata[h] = Homolog.get(*attrs)

        return metadata

    def homolog(self, id):
        """ Return Homolog with id. """
        return self.homologs[id]

    def map(self, attr1, attr2=None):
        """ Return mapping between two attributes in homolog set, OR
            if no second attribute is give, map attr1 to whole homolog.

            attr2 can be a list of other attributes -- which will return
                a list of those attributes mapped to attr1.
        """
        m = dict()

        # If no second attribute is give, mapping is between first attribute
        # and the homolog object
        if attr2 is None:
            for id, h in self._homologs.items():
                m[getattr(h, attr1)] = h

        # else, mapping from one attribute to another
        else:

            # If attr2 is a list, return a list in mapping
            if isinstance(attr2,list):
                for id, h in self._homologs.items():
                    m[getattr(h, attr1)] = [getattr(h, a) for a in attr2]

            # Else just return a single attr.
            else:
                for id, h in self._homologs.items():
                    m[getattr(h, attr1)] = getattr(h, attr2)

        return m

    def subset(self, ids, inplace=False):
        """Subset the HomologSet. By default, it reduces
        the HomologSet size in place. If `inplace`=False,
        subset HomologSet is returned as new_object
        """
        if inplace:
            homologs_to_rm = [h for h in self.homologs if h not in ids]
            self.rm(homologs_to_rm)
        else:
            homologs = []
            for i in ids:
                homologs.append(self.homologs[i])

            return HomologSet(homologs=homologs)

    def renumber(self):
        """ Renumber the ID numbers for homologs in the set, starting at
            0 to len(homologs).
        """
        n = len(self.homologs)
        for i in range(n):
            unique_id = "XX%08d" % i
            self._homologs[i].id = unique_id

    def add(self, homologs={}, **kw_homologs):
        """ Append a list of homolog objects to the set.

            NOTE: does not renumber the homolog set. this must be called
            manually.
        """
        # If a single homolog is given, format it into a list
        # for loop below.
        if isinstance(homologs, list) == False:
            homologs = [homologs]

        # Set homolog as attribute of the HomologSet object
        for h in homologs:
            setattr(self, h.id, h)

        # Add homolog to homologs
        for h in homologs:
            self._homologs[h.id] = h

    def rm(self, ids):
        """ Remove a list of homologs from set of homologs.

            NOTE: this does not renumber the homolog set. must be done manually.

            Also, note that this changes the homolog set in place!
        """
        # If a single id is given, format it into a list for loop below.
        if isinstance(ids,list) == False:
            ids = [ids]

        for i in ids:
            # Remove from homologs dict
            del self._homologs[i]

            # Remove Homolog attribute from object
            delattr(self, i)

    def cluster(self,
        redund_cutoff=0.99,
        tmp_file_suffix="oB_cdhit",
        word_size=5,
        cores=1,
        keep_tmp=False,
        accession=(),
        positive=(),
        negative=("putative","hypothetical", "unnamed", "possible", "predicted",
                    "unknown", "uncharacterized","mutant", "isoform"),
        inplace=False
        ):
        """ Reduce sequence redundancy in the HomologSet using CDHIT.

        By default, this method returns a new HomologSet

        """
        # Write out the fasta file with a unique name for each sequence that goes
        # >0, >1...>N.  Those numbers point to the index of the sequence in
        # homolog_list.

        # Don't do anything for empty list
        if len(self.homologs) == 0:
            print("Warning: empty list passed to cdhit!  Ignoring.")
            return self

        # Ranks homologs.
        try:
            rank_homologs(self, accession=accession, positive=positive, negative=negative)
        except:
            print("""WARNING: No defline in Homologs?""")

        # Create a temporary fasta file from homologs as input to CDHIT.
        fname_prefix = "cluster-%s" % str(redund_cutoff)
        fname = fname_prefix+".fasta"
        self.Write.fasta(fname=fname, tags=["id"])

        cdhit.run(fname_prefix,
            redund_cutoff=redund_cutoff,
            tmp_file_suffix=tmp_file_suffix,
            word_size=word_size,
            cores=cores,
            keep_tmp=keep_tmp,
            accession=accession,
            positive=positive,
            negative=negative
        )
        print("""Clustering finished.""")

        # Now parse the output of cdhit and grab members of clusters with the
        # lowest rank
        f = open("%s_cdhit.clstr" % fname_prefix,'r')

        # Get a id-to-rank mapping dict from homolog_set
        id_rank = self.map("id", "rank")

        subset_ids = []
        in_cluster = []
        line = f.readline()
        while line != "":

            # If we are starting a new cluster
            if line.startswith(">"):

                # ... and this is not the first cluster
                if in_cluster != []:

                    # Take the member of in_cluster with the minimum rank
                    ranks = [id_rank[homolog_id] for homolog_id in in_cluster]
                    best_id = in_cluster[ranks.index(min(ranks))]
                    subset_ids.append(best_id)

                in_cluster = []

            # If this is not a blank line, record the seq_id in in_cluster
            elif line[0] in string.digits:
                seq_id = line.split(">")[1][:10]
                in_cluster.append(seq_id)

            # Read the next line
            line = f.readline()

        # Grab the last cluster
        ranks = [id_rank[homolog_id] for homolog_id in in_cluster]
        best_id = in_cluster[ranks.index(min(ranks))]
        subset_ids.append(best_id)
        f.close()

        # Delete temporary files
        if not keep_tmp:
            try:
                os.remove("%s.fasta" % fname_prefix)
                os.remove("%s_cdhit" % fname_prefix)
                os.remove("%s_cdhit.clstr" % fname_prefix)
                os.remove("%s_cdhit.bak.clstr" % fname_prefix)
            except:
                pass

        # Subset the HomologSet with ids pulled from clusters.
        if inplace:
            self.subset(subset_ids, inplace=inplace)
            print("""%s sequences in HomologSet""" % len(self.homologs))
        else:
            hs = self.subset(subset_ids, inplace=inplace)
            return hs
