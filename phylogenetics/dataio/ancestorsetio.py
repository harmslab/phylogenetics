import base
import dendropy

from .formats import rst

class Write(base.Write):

    def __init__(self, AncestorSet):
        """ """
        self._AncestorSet = AncestorSet

    def _ancestor_set_to_ancestor_data(self):
        """ """

    def fasta(self):
        """ Return AncestorSet as fasta."""

    def csv(self):
        """ """
        pass


class Read(base.Read):

    def __init__(self, AncestorSet):
        """ Read AncestorSet data from file and add to Ancestor"""
        self._AncestorSet = AncestorSet

    def _ancestor_data_to_ancestor_set(self, tree, ancestor_data):
        """ Take the ancestor tree, and ancestor data dictionary, and map it
            to the tree attached to AncestorSet object.

            Must take a newick tree format.

            Ancestor data:
            -------------
            {"ZZ00000000":
                {0:
                    "A":0.0, "C":0.0, ...
                },
                {1:
                    "A":0.8, "C":0.0, ...
                }
                ...,
             "ZZ00000001":
                {0:
                    "A":0.1, "C":0.6, ...
                },
                {1:
                    "A":0.3, "C":0.0, ...
                }
                ...,
            ...
            }
        """
        # Read newick file into dendropy tree object
        new_tree = dendropy.datamodel.treemodel.Tree.get_from_string(tree, schema="newick")
        old_tree = self._AncestorSet._Tree._DendroPyTree

        new_nodes = new_tree.internal_nodes()
        old_nodes = old_tree.internal_nodes()

        # Sanity check that nodes in loaded tree match old nodes
        if len(new_nodes) != len(old_nodes):
            raise Exception("""Nodes in loaded tree do not match tree topology of\
            of AncestorSet tree.
            """)

        # Iterate through internal nodes and bind reconstruction data.
        for i in range(len(new_nodes)):
            # Get ancestor's label from reconstructed tree
            anc_key = int(new_nodes[i].label)
            # Get reconstruction data from tree
            anc_data = ancestor_data[anc_key]

            # Bind reconstructed ancestor to AncestorSet
            Ancestor = getattr(self._AncestorSet, old_nodes[i].label)
            Ancestor.sites = anc_data

        return self._AncestorSet


    @base.read_from_file
    def rst(self, data):
        """ Read ancestor set from paml file. """
        tree, ancestor_data = rst.read(data)
        self._ancestor_data_to_ancestor_set(tree, ancestor_data)
        return self._AncestorSet
