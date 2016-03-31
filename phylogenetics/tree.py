import dendropy

class Tree(dendropy.datamodel.treemodel.Tree):
    """
        Subclass of Dendropy's Tree object including extra methods for
        unifying with homologsets.
    """

    @classmethod
    def get_from_homologset(cls, homologset):
        """
            Get from src.
        """
        # Use the get method from Super class
        instance = super(Tree, cls).get(data=homologset.tree, schema="newick")

        # Create an instance
        instance._homologset = homologset

        # Add metadata to tree homologs and nodes
        instance._tip_metadata()
        instance._node_metadata()
        instance._edge_metadata()

        # Return instance
        return instance

    def prune_homologs(self, ancestor):
        """
            Prune the subtree below named ancestral node.

            Argument:
            --------
            ancestor: str
                Ancestor label
        """
        # If given a number for the ancestor, convert to ancestor label.
        if type(ancestor) == int:
            ancestor = "Anc" + str(ancestor)

        # Find the node with the given label
        node = self.find_node(lambda n: n.taxon.label==ancestor)

        if node is None:
            raise Exception(""" No ancestor with the given name. """)

        # Get list of child homologs for this node
        homologs = node.leaf_nodes()

        # Read their ids
        ids = [h.taxon.label for h in homologs]

        # Remove homologs from homolog_set
        self._homologset.rm_homologs(ids)

    def _tip_metadata(self, attributes=("accession", "organism", "alignedlen")):
        """
            Set the metadata of the tips of the tree.
        """
        # Iterate through the attributes given
        for a in attributes:
            # If the attribute doesn't exist, skip it
            try:
                mapping = self._homologset.get_map("id", a)

                # Add that attribute to tree's taxon metadata
                for taxon in self.taxon_namespace:
                    id = taxon.label
                    taxon.annotations.add_new(a, mapping[id])
            except:
                pass


    def _node_metadata(self):
        """
            Get metadata for internal nodes
        """
        # Get the nodess
        nodes = self.nodes()

        # Iterate through nodes and look for nodes with node labels
        for i in range(len(nodes)):

            # If the node taxon is None, then create Taxon object
            if nodes[i].taxon is None:
                # Add taxon object
                nodes[i].annotations.add_new("name","Anc" + str(i))

    def _edge_metadata(self):
        """
            Set the metadata of the edges
        """
        # Get the nodess
        edges = self.edges()

        # Iterate through nodes and look for nodes with node labels
        for i in range(len(edges)):
            edges[i].annotations.add_new("name", str(i))


    def color_homologs(self, ids, color):
        """
            Color (hexadecimal code only) nodes with given ids.
        """
        # Find the nodes given
        nodes = self.find_node(lambda n: n.taxon.label in ids)

        # Add color as metadata
        for n in nodes:
            n.taxon.annotations.add_new("!color", color)


def add_tree_to_homologset(homologset, newick):
    """
        Adds a Tree object to homolog set.
    """
    homologset.Tree = Tree.get_from_homologset(homologset)
    return homologset


def ancestral_labels(dendrotree):
    """
        Labels the ancestral nodes with 'Anc##'.

        Does this by looking for nodes with no taxa. Not the most
        intelligent way to do this yet.

        Arguments:
        ---------
        dendrotree: dendropy.datamodel.treemodel
            Dendropy Tree object

        Returns:
        -------
        dendrotree:
            Tree with new taxon labels.
    """
    # Get the nodess
    nodes = dendrotree.nodes()

    # Iterate through nodes and look for nodes with node labels
    for i in range(len(nodes)):

        # If the node taxon is None, then create Taxon object
        if nodes[i].taxon is None:
            # Add taxon object
            nodes[i].annotations.add_new("ancestor_label","Anc" + str(i))

    return dendrotree
