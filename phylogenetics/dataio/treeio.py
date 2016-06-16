import .base
class Write(base.Write):
    """Write the tree out."""
    def __init__(self, Tree):
        self._Tree = Tree

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
        output = self._Tree.Dendropy.as_string(schema="nexus",
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
        output = self._Tree.Dendropy.as_string(schema="newick",
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

    def nexus(self):
        """ """
