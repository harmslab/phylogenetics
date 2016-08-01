# Tree object for HomologSets
from __future__ import absolute_import

import dendropy
import copy
from phylogenetics.dataio import treeio

class Tree(object):
    """Subclass of Dendropy's Tree object including extra methods for
    connectiong to HomologSet object.


    """
    def __init__(self, HomologSet, tree="", stats={}):
        # Stats about tree construction in dictionary.
        self.stats = stats
        self._HomologSet = HomologSet

        # Add Write module to tree
        self.Write = treeio.Write(self)
        self.Read = treeio.Read(self)

    @property
    def DendroPy(self):
        """ """
        try:
            return self._DendroPy
        except AttributeError:
            raise Exception("""\
            DendroPy object was never initialized. Try reading in a tree via the\
            `Read` subobject attached this instance.\
            """)

    def _construct(self):
        """Construct a Dendropy object
        """
        # Bind nodes to homologs
        self._nodes_to_homologs()
        self._label_nodes()
        self._label_edges()

    def _nodes_to_homologs(self):
        """ Point nodes to homolog object and vice versa."""
        # Iterate through nodes in tree.
        for node in self.DendroPy.leaf_nodes():

            # Find nodes that represent tips of the tree.
            #if node.taxon is not None:
            # Get the id of the homolog
            id = node.taxon.label

            # Get Homolog object for this node
            Homolog = getattr(self._HomologSet, id)

            # Bind Homolog to node
            node.Homolog = Homolog

            # Bind node object to Homolog
            setattr(Homolog, "node", node)


    def label_tips(self, *attrs, clear=False):
        """Add Homolog Attributes as tree tip labels. Pass then name of the
        attributes as arguments into this method.

        Set clear to True if you'd like to clear any previously added labels
        """
        for node in self.DendroPy.leaf_nodes():
            # Clear old annotations?
            if clear:
                node.annotations.clear()
            for a in attrs:
                node.annotations.add_new(a, getattr(node.Homolog, a))

    def _label_nodes(self):
        """ Label the internal nodes of the tree (a.k.a. ancestors). """
        i = 0
        for node in self.DendroPy.internal_nodes():
            # Get the labels given by tree software.
            label = node.label
            node.score = label
            # Define a unique label for internal label
            node.label = "ZZ%08d" % i
            i += 1

    def _label_edges(self):
        """Label the edges of the tree.
        """
        i = 0
        for edge in self.DendroPy.edges():
            # Give the edge a unique label
            edge.label = "YY%08d" % i
            i += 1

    def subtree(self, node_name):
        """ Get a subtree. """
        # Find the node with the given label
        node = self.DendroPy.find_node(lambda n: n.annotations.get_value("name")==node_name)

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
        homologset = self._HomologSet.subset(ids)
        # 5. Add the trees to homologetset.
        homologset.tree = tree
        homologset.Tree = Tree

        return homologset

    def reroot(self, item_name):
        """Reroot a tree on a given item.

        Parameters
        ----------
        item_name : str
            Name of object in tree to reroot the tree. Can be either a ancestor,
            edge, or taxa.
        """
        if type(item_name) == str:
            # Item dictionary
            methods = {
                "XX" : (self.DendroPy.nodes, self.DendroPy.reroot_at_node),
                "YY" : (self.DendroPy.edges, self.DendroPy.reroot_at_edge),
                "ZZ" : (self.DendroPy.internal_nodes, self.DendroPy.reroot_at_node)
            }
            # Retrieve method for finding reroot object
            search_method = methods[item_name[0:2]][0]
            reroot_method = methods[item_name[0:2]][1]

        # Iterate through all items
        for thing in search_method():
            # Find the item with the given label
            if thing.label == item_name:
                item = thing

        # Reroot tree at item
        reroot_method(item)

    def prune(self, id):
        """Prune Node in Tree object. Also removes Homolog from HomologSet.
        """
        # Get homolog object
        Homolog = getattr(self._HomologSet, id)

        # Remove node (and children) from tree.
        self.DendroPy.prune_subtree(Homolog.node)

        # Remove Homolog from HomologSet
        self._HomologSet.rm(id)
