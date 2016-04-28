# In the future, this will be the top level object that contains all Subobjects
# important for doing a phylogenetics project

class PhylogeneticProject(object):
    """Container object for managing all data for a phylogenetics project.
    """
    def __init__(self, **kwargs):
        """
            Base Object for managing phylogenetic projects
        """
        pass

    def add(self, item):
        """Add data to project.
        """
        # Get type of item
        item_type = type(item)

        # possible objects to add
        items = {
            str(type(HomologSet)): self._add_HomologSet,
            str(type(Alignment)): self._add_Alignment,
            str(type(Tree)): self._add_Tree,
            str(type(Reconstruction)): self._add_Reconstruction
        }

        # Find item type in set of possible items
        adding_method = items[item_type]

        # Add that item to project
        adding_method(item)

    @classmethod
    def load(self, fname):
        """Load a phylogenetics project from file.
        """

    def _add_HomologSet(self, HomologSet):
        """Add HomologSet Set to PhylogeneticsProject object."""
        # Set the HomologSet object
        self.HomologSet = HomologSet
        # Expose the align method of this object to user
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

    def _align(self, fname="alignment.fasta", rm_tmp=True, quiet=False):
        """ Multiple sequence alignment of the HomologSet.

            Currently, only option is to use MSAProbs.
        """
        # Write out alignment file
        self.HomologSet.Write.fasta(fname="alignment.fasta")

        # Run the alignment with MSAProbs
        output_fname = msaprobs.run(fasta_fname="alignment", rm_tmp=rm_tmp)

        # Attach an alignment object to HomologSet
        self._add_Alignment(Alignment(self))

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
        self._add_Tree = Tree(self, tree, stats=stats)


    def _reconstruct(self):
        """ Resurrect Ancestors on Tree.
        """
        # Bind Ancestor Objects to each internal node.
        ancestors = []
        for node in self.Tree._DendroPyTree.internal_nodes():
            id = node.label
            ancestors.append( Ancestor(id, self.Tree))

        # Bind an AncestorSet object to HomologSet
        self.AncestorSet = AncestorSet(self.Tree, ancestors=ancestors)
        self.AncestorSet._nodes_to_ancestor()


        



















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
