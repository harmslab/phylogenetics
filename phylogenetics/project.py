# In the future, this will be the top level object that contains all Subobjects
# important for doing a phylogenetics project
import os

# import objects to bind to Project class
from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.alignment import Alignment
from phylogenetics.tree import Tree
from phylogenetics.ancestors import Ancestor, AncestorSet
from phylogenetics.reconstruction import Reconstruction

from phylogenetics.dataio import projectio

# imports for running external tools.
from .exttools import (cdhit,
                        msaprobs,
                        phyml,
                        paml)


class Project(object):
    """Container object for managing all data for a phylogenetics project.

    Optional arguments can be passed into the Project class. These must be
    phylogenetic objects from this package (i.e. HomologSet, Alignment, Tree,
    Reconstruction, AncestorSet, etc.)
    """
    def __init__(self, *args):
        # Bind Reading module to class
        self.Read = projectio.Read(self)
        # Add any objects that were given to Project
        for a in args:
            self.add(a)

    def add(self, item):
        """Add data to project.
        """
        # possible objects to add
        items = {
            HomologSet: self._add_HomologSet,
            Alignment: self._add_Alignment,
            Tree: self._add_Tree,
            Reconstruction: self._add_Reconstruction,
            AncestorSet: self._add_AncestorSet
        }

        # Find item type in set of possible items
        adding_method = items[item.__class__]

        # Add that item to project
        adding_method(item)

    def download(self, ids, email):
        """ Download a set of Homologs"""
        hs = HomologSet()
        self._add_HomologSet( hs )
        self.HomologSet.download(ids, email)

    def _add_HomologSet(self, HomologSet):
        """Add HomologSet Set to PhylogeneticsProject object."""
        # Set the HomologSet object
        self.HomologSet = HomologSet
        # Expose the align method of this object to user
        setattr(self, "cluster", self._cluster)
        setattr(self, "align", self._align)

    def _add_Alignment(self, Alignment):
        """Add Alignment to PhylogeneticsProject object."""
        # Set the Alignment object
        self.Alignment = Alignment
        # Expose the tree methods of this project object
        setattr(self, "tree", self._tree)

    def _add_Tree(self, Tree):
        """Add Tree to PhylogeneticsProject object."""
        # Set the Tree object of project
        self.Tree = Tree
        # Expose the reconstruction methods of this project object
        setattr(self, "reconstruct", self._reconstruct)

    def _add_Reconstruction(self, Reconstruction):
        """Add Reconstruction to PhylogeneticsProject object."""
        self.Reconstruction = Reconstruction

    def _add_AncestorSet(self, AncestorSet):
        """Add a AncestorSet object to PhylogeneticsProject object."""
        self.AncestorSet = AncestorSet

    def _cluster(self,
        redund_cutoff=0.99,
        tmp_file_suffix="oB_cdhit",
        word_size=5,
        cores=1,
        keep_tmp=False,
        accession=(),
        positive=(),
        negative=("putative","hypothetical", "unnamed", "possible", "predicted",
                    "unknown", "uncharacterized","mutant", "isoform"),
        inplace=True
        ):
        """Remove redundant sequences from HomologSet.
        """
        self.HomologSet.cluster(
            redund_cutoff=redund_cutoff,
            tmp_file_suffix=tmp_file_suffix,
            word_size=word_size,
            cores=cores,
            keep_tmp=keep_tmp,
            accession=accession,
            positive=positive,
            negative=negative
        )


    def _align(self, fname="alignment.fasta", rm_tmp=True, quiet=False):
        """ Multiple sequence alignment of the HomologSet.

            Currently, only option is to use MSAProbs.
        """
        # Write out alignment file
        self.HomologSet.Write.fasta(fname="alignment.fasta")

        # Run the alignment with MSAProbs
        output_fname = msaprobs.run(fasta_fname="alignment", rm_tmp=rm_tmp)

        # Attach an alignment object to HomologSet
        self._add_Alignment(Alignment(self.HomologSet))

        # Read alignment from output fasta and manage with Alignment object
        self.Alignment.Read.fasta(fname=output_fname)

        # Let us know when finished
        if quiet is False:
            print("Alignment finished.")

        # Remove fasta file.
        if rm_tmp:
            os.remove(output_fname)

    def _tree(self, **kwargs):
        """ Compute the maximum likelihood phylogenetic tree from
            aligned dataset.

        """
        # Write the HomologSet out as a phylip.
        self.Alignment.Write.phylip(fname="ml-tree.phy")

        # Run phyml and parse results.
        tree, stats = phyml.run("ml-tree", **kwargs)

        # Add Tree object to HomologSet
        self._add_Tree(Tree(self.HomologSet, tree, stats=stats))


    def _reconstruct(self):
        """ Resurrect Ancestors on Tree.
        """
        # Bind Ancestor Objects to each internal node.
        ancestors = []
        for node in self.Tree._DendroPyTree.internal_nodes():
            id = node.label
            ancestors.append( Ancestor(id, self.Tree))

        # Bind an AncestorSet object to HomologSet
        self._add_AncestorSet( AncestorSet(self.Tree, ancestors=ancestors) )
        self.AncestorSet._nodes_to_ancestor()

        seqfile = "asr-alignment.fasta"
        outfile = "asr-output"
        treefile = "asr-tree.nwk"

        # Prepare input files for PAML
        self.Tree._DendroPyTree.write(path=treefile, schema="newick", suppress_internal_node_labels=True)
        self.Alignment.Write.fasta(fname=seqfile)

        # Construct a paml job
        paml_job = paml.CodeML(
            seqfile=seqfile,
            outfile=outfile,
            treefile=treefile,
            fix_alpha=True,
            alpha=self.Tree.stats["Gamma shape parameter"],
        )

        reconstruction = Reconstruction(self.Alignment, self.Tree, self.AncestorSet, paml_job)
        self._add_Reconstruction( reconstruction )

        # Run the PAML job
        self.Reconstruction.paml_job.run()

        # Read the paml output and bind data to tree
        self.AncestorSet.Read.rst(fname="rst")

        # Infer gaps.
        self.Reconstruction.infer_gaps()
