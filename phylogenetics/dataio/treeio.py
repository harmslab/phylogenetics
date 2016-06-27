from . import base
from dendropy.datamodel.treemodel import Tree as DendropyTree

class Write(base.Write):
    """Write the tree out."""
    def __init__(self, Tree):
        self._Tree = Tree

    def _object_to_data(self):
        """Write object into metadata"""
        data = self.newick(
            suppress_leaf_taxon_labels=False,
            suppress_leaf_node_labels=True,
            suppress_internal_taxon_labels=False,
            suppress_internal_node_labels=False,
            suppress_rooting=False,
            suppress_edge_lengths=False,
            unquoted_underscores=False,
            preserve_spaces=False,
            store_tree_weights=False,
            taxon_token_map=None,
            suppress_annotations=True,
            annotations_as_nhx=False,
            suppress_item_comments=True,
            node_label_element_separator=' ',
            node_label_compose_fn=None,
            edge_label_compose_fn=None,
            real_value_format_specifier='',
            ignore_unrecognized_keyword_arguments=False,
        )
        return data

    @base.write_to_file
    def nexus(self,
        suppress_leaf_taxon_labels=False,
        suppress_leaf_node_labels=True,
        suppress_internal_taxon_labels=False,
        suppress_internal_node_labels=False,
        suppress_rooting=False,
        suppress_edge_lengths=False,
        unquoted_underscores=False,
        preserve_spaces=False,
        store_tree_weights=False,
        taxon_token_map=None,
        suppress_annotations=False,
        annotations_as_nhx=False,
        suppress_item_comments=True,
        node_label_element_separator=' ',
        node_label_compose_fn=None,
        edge_label_compose_fn=None,
        real_value_format_specifier='',
        ignore_unrecognized_keyword_arguments=False,):
        """Write tree to nexus format -- using DendroPy machinery.

        For argument options, see DendroPy documentation.
        """
        output = self._Tree.DendroPy.as_string(schema="nexus",
            suppress_leaf_taxon_labels=False,
            suppress_leaf_node_labels=True,
            suppress_internal_taxon_labels=False,
            suppress_internal_node_labels=False,
            suppress_rooting=False,
            suppress_edge_lengths=False,
            unquoted_underscores=False,
            preserve_spaces=False,
            store_tree_weights=False,
            taxon_token_map=None,
            suppress_annotations=False,
            annotations_as_nhx=False,
            suppress_item_comments=True,
            node_label_element_separator=' ',
            node_label_compose_fn=None,
            edge_label_compose_fn=None,
            real_value_format_specifier='',
            ignore_unrecognized_keyword_arguments=False
        )
        return output


    @base.write_to_file
    def newick(self,
        suppress_leaf_taxon_labels=False,
        suppress_leaf_node_labels=True,
        suppress_internal_taxon_labels=False,
        suppress_internal_node_labels=False,
        suppress_rooting=False,
        suppress_edge_lengths=False,
        unquoted_underscores=False,
        preserve_spaces=False,
        store_tree_weights=False,
        taxon_token_map=None,
        suppress_annotations=False,
        annotations_as_nhx=False,
        suppress_item_comments=True,
        node_label_element_separator=' ',
        node_label_compose_fn=None,
        edge_label_compose_fn=None,
        real_value_format_specifier='',
        ignore_unrecognized_keyword_arguments=False,):
        """Write tree to newick format -- using DendroPy machinery.

        For argument options, see DendroPy documentation.
        """
        output = self._Tree.DendroPy.as_string(schema="newick",
            suppress_leaf_taxon_labels=False,
            suppress_leaf_node_labels=True,
            suppress_internal_taxon_labels=False,
            suppress_internal_node_labels=False,
            suppress_rooting=False,
            suppress_edge_lengths=False,
            unquoted_underscores=False,
            preserve_spaces=False,
            store_tree_weights=False,
            taxon_token_map=None,
            suppress_annotations=True,
            annotations_as_nhx=False,
            suppress_item_comments=True,
            node_label_element_separator=' ',
            node_label_compose_fn=None,
            edge_label_compose_fn=None,
            real_value_format_specifier='',
            ignore_unrecognized_keyword_arguments=False
        )
        return output

class Read(base.Read):
    """"""
    def __init__(self, Tree):
        self._Tree = Tree

    def _data_to_object(self, data, schema=None):
        """Attaches a DendroPy tree object"""
        schemas = [
            "newick",
            "nexus"
        ]
        # If schemas is not given, try each schema until one works.
        if schema is None:
            for s in schemas:
                # Try different schemas until a tree works
                try:
                    tree = DendropyTree.get(data=data, schema=s)
                except:
                    pass
        else:
            tree = DendropyTree.get(data=data, schema=schema)

        # Check that a tree was made.
        try:
            self._Tree._DendroPy = tree
            self._Tree._construct()
        except NameError:
            raise Exception("""Tree data doesn't seem to be in a format that DendroPy can read.""")

        return self._Tree

    def newick(self, data):
        """Read a tree from nexus string."""
        return self._data_to_object(data=data, schema="newick")

    def nexus(self, data):
        Tree = DendropyTree.get(data=data, schema="nexus")
        return self._data_to_object(Tree)
