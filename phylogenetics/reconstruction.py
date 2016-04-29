import dendropy

class Reconstruction(object):
    """Object for doing ancestral sequence reconstruction
    """
    def __init__(self, Alignment, Tree, AncestorSet, paml_job, gaps_infered=True):
        """ """
        self._Alignment = Alignment
        self._Tree = Tree
        self._AncestorSet = AncestorSet

        self.gaps_infered = gaps_infered
        self.paml_job = paml_job

    def infer_gaps(self):
        """Estimate gaps for ancestral nodes.


        Handling Gaps
        -------------
        Uses Fitch Parsimony to determine sites in ancestral sequences that
        are likely gaps. The posterior probability of these sites are ignored
        when calculating the average posterior probability of the ancestor.
        """
        taxa = self._Tree._DendroPyTree.taxon_namespace

        # Build a Sequence data matrix from Dendropy
        data = dendropy.ProteinCharacterMatrix.get(
            data=self._Alignment.Write.fasta(),
            schema="fasta",
            taxon_namespace=taxa
        )

        tree = self._Tree._DendroPyTree

        # Get the alphabet of Dendropy's ProteinCharacterMatrix
        alphabet = data.state_alphabets[0].symbols

        # Construct a map object between sequence data and tree data.
        taxon_state_sets_map = data.taxon_state_sets_map(gaps_as_missing=False)

        # Fitch algorithm to determine placement of gaps
        dendropy.model.parsimony.fitch_down_pass(tree.postorder_node_iter(),
                taxon_state_sets_map=taxon_state_sets_map)
        dendropy.model.parsimony.fitch_up_pass(tree.preorder_node_iter())

        # Iterate through each ancestor and find their gaps
        for id, Ancestor in self._AncestorSet.ancestors.items():
            gap_indices = []

            # Iterate through the max parsimony at each set for the ancestor to
            # find any gaps that should be placed
            for i in range(len(Ancestor.node.state_sets)):
                site_list = list(Ancestor.node.state_sets[i])

                # Find gaps in list and add to gap array
                if len(site_list) == 1 and alphabet[site_list[0]] == "-":
                    gap_indices.append(i)

            # Bind gap indices to Ancestor
            if gap_indices != []:
                Ancestor._set_gaps(gap_indices)
