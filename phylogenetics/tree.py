import dendropy
import copy

from phylogenetics.homologs import HomologSet

class Ancestor:

    def __init__(self, Tree, data):

        ## Quality control name!
        self._Tree = Tree
        self._data = data

    @property
    def posteriors(self):
        return self._posteriors

    @property
    def mlsequence(self):
        return self._mlsequence


class AncestorSet(object):

    def __init__(self, Tree, Ancestors=None):
        """ Object that holds a set of ancestor object. """

        self._Tree = Tree

        if Ancestors is not None:
            for a in Ancestors:
                self.add(a)

    def add(self, Ancestor):
        """ Add Ancestor to set. """
        setattr(self, Ancestors.name, Ancestors)

    def rm(self, Ancestor):
        """ Remove the Ancestor from a set. """
        delattr(self, Ancestor.name)

class Tree(dendropy.datamodel.treemodel.Tree):
    """
        Subclass of Dendropy's Tree object including extra methods for
        connectiong to HomologSet object.
    """

    @classmethod
    def get_from_homologset(cls, homologset, attributes=("accession")):
        """
            Get from src.
        """
        # Use the get method from Super class
        instance = super(Tree, cls).get(data=homologset.tree, schema="newick")

        # Create an instance
        instance._homologset = homologset

        # Add metadata to tree homologs and nodes
        instance.add_tip_metadata(attributes)
        instance.add_node_metadata()
        instance.add_edge_metadata()

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
        node = self.find_node(lambda n: n.annotations.get_value("name")==ancestor)

        if node is None:
            raise Exception(""" No ancestor with the given name. """)

        # Get list of child homologs for this node
        homologs = node.leaf_nodes()

        # Read their ids
        ids = [h.taxon.label for h in homologs]

        # Remove homologs from homolog_set
        self._homologset.rm_homologs(ids)

        # Prune subtree
        self.prune_subtree(node)

    def subtree(self, node_name):
        """ Get a subtree. """
        # Find the node with the given label
        node = self.find_node(lambda n: n.annotations.get_value("name")==node_name)

        if node is None:
            raise Exception(""" No ancestor with the given name. """)

        # Get list of child homologs for this node
        tips = node.leaf_nodes()

        # Read their ids
        ids = [h.taxon.label for h in tips]

        non_interesting_ids = self._homologset.id_list

        # Remove interesting ids
        for id in ids:
            non_interesting_ids.remove(id)

        # Sadly, this logic is a little confusing and inefficient --
        # not sure how to make it better.
        # 1. Clone tree
        Tree = copy.deepcopy(self)
        # 2. Prune tree with labels
        Tree.prune_taxa_with_labels(non_interesting_ids)
        # 3. Get newick string
        tree = Tree.as_string(schema="newick")
        # 4. Construct a new homologset
        homologset = self._homologset.subset(ids)
        # 5. Add the trees to homologetset.
        homologset.tree = tree
        homologset.Tree = Tree

        return homologset


    def add_tip_metadata(self, attributes):
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


    def add_node_metadata(self):
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
            else:
                # Initialize annotations
                nodes[i].annotations.add_new("name","Tip" + str(i))

    def add_edge_metadata(self):
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
