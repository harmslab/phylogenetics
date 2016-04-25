# API for working with Homolog sets in a phylogenetics project
import os

# Reading and Writing modules to HomologSet
import phylogenetics.dataio.homologio as homologio
import phylogenetics.dataio.homologsetio as homologsetio

# Subobjects bound to HomologSet
from phylogenetics.alignment import Alignment
from phylogenetics.tree import Tree
from phylogenetics.ancestors import Ancestor, AncestorSet

# Import external tools
from phylogenetics.exttools import (cdhit,
                                msaprobs,
                                phyml,
                                paml)


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

    @property
    def attrs(self):
        """ """
        return self._attrs

    @property
    def seqlen(self):
        """ Get the length of sequence """
        return len(self.sequence)

    def add_alignment(self, alignment):
        """ Add alignment to homolog """
        # Check if this is the first alignment
        if hasattr(self, "latest_align") is False:
            # If first alignment, add new attribute to homolog.
            self.addattr("latest_align", alignment)

        # Else move old alignment to new attribute and update latest alignment.
        else:
            # Look for latest alignment
            counter = 0
            attr = "align0"
            while hasattr(self, attr) is True:
                attr = "align%d" % counter
                counter += 1

            # move old alignment to new attribute.
            self.addattr(attr, self.latest_align)

            # add new alignment to latest_align attribute
            self.latest_align = alignment

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

    @classmethod
    def load(cls, fname):
        """ Read a HomologSet. """
        with open(fname, "rb") as f:
            homologset = pickle.load(f)
        return homologset

    @property
    def homologs(self):
        """ Get homolog set. """
        return self._homologs

    @property
    def list_ids(self):
        """ Return ID list. """
        return [h.id for h in self._homologs]

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

    def subset(self, ids):
        """ Get a subset of the homologs from an ID list. """
        mapping = self.HomologSetmap("id")

        homologs = []
        for id in ids:
            homologs.append(mapping[id])

        return HomologSet(homologs=homologs)

    def renumber(self):
        """ Renumber the ID numbers for homologs in the set, starting at
            0 to len(homologs).
        """
        n = len(self.homologs)
        for i in range(n):
            unique_id = "XX%08d" % i
            self._homologs[i].id = unique_id

    def add(self, homologs):
        """ Append a list of homolog objects to the set.

            NOTE: does not renumber the homolog set. this must be called
            manually.
        """
        # If a single homolog is given, format it into a list
        # for loop below.
        if isinstance(homologs,list) == False:
            homologs = [homologs]

        # Set homolog as attribute of the HomologSet object
        for h in homologs:
            setattr(self, h.id, h)

        # Add homolog to homologs
        for h in homologs:
            self._homologs[h.id] = h


    def cluster(self):
        """ Reduce any redundancy in the HomologSet using CDHIT.

        """
        pass


    def rm(self, ids):
        """ Remove a list of homologs from set of homologs.

            NOTE: this does not renumber the homolog set. must be done manually.

            Also, note that this changes the homolog set in place!
        """
        # If a single id is given, format it into a list for loop below.
        if isinstance(ids,list) == False:
            ids = [ids]

        for id in ids:
            # Remove from homologs dict
            del self._homologs[id]

            # Remove Homolog attribute from object
            delattr(self, id)

    def align(self, fname="alignment.fasta", rm_tmp=True, quiet=False):
        """ Multiple sequence alignment of the HomologSet.

            Currently, only option is to use MSAProbs.
        """
        # Write out alignment file
        self.Write.fasta(fname="alignment.fasta")

        # Run the alignment with MSAProbs
        output_fname = msaprobs.run(fasta_fname="alignment", rm_tmp=rm_tmp)

        # Attach an alignment object to HomologSet
        self.Alignment = Alignment(self)

        # Read alignment from output fasta and manage with Alignment object
        self.Alignment.Read.fasta(fname=output_fname)

        # Let us know when finished
        if quiet is False:
            print("Alignment finished.")

        # Remove fasta file.
        if rm_tmp:
            os.remove("%s.fasta" % output_fname)

    def tree(self, **kwargs):
        """ Compute the maximum likelihood phylogenetic tree from
            aligned dataset.

        """
        # Write the HomologSet out as a phylip.
        self.Alignment.Write.phylip(fname="ml-tree.phy")

        # Run phyml and parse results.
        tree, stats = phyml.run("ml-tree", **kwargs)

        # Add Tree object to HomologSet
        self.Tree = Tree(self, tree, stats=stats)


    def reconstruct(self):
        """ Resurrect Ancestors on Tree.
        """
        # Bind Ancestor Objects to each internal node.
        ancestors = []
        for node in self.Tree._DendroPyTree.internal_nodes():
            id = node.label
            ancestors.append( Ancestor(id, self.Tree))

        # Bind an AncestorSet object to HomologSet
        self.AncestorSet = AncestorSet(self.Tree, ancestors=ancestors)

        seqfile = "asr-alignment.fasta"
        outfile = "asr-output"
        treefile = "asr-tree.nwk"

        # Prepare input files for PAML
        self.Tree._DendroPyTree.write(path=treefile, schema="newick", suppress_internal_node_labels=True)
        self.Alignment.Write.fasta(fname=seqfile)

        # Construct a paml job
        pamljob = paml.CodeML(
            seqfile=seqfile,
            outfile=outfile,
            treefile=treefile,
            fix_alpha=True,
            alpha=self.Tree.stats["Gamma shape parameter"],
        )

        # Run the PAML job
        pamljob.run()

        # Read the paml output and bind data to tree
        self.AncestorSet.Read.rst(fname="rst")
