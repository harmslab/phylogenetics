import dendropy
from . import handlers

class DendroPyError(Exception):
    """DendroPy object was never initialized. Try reading in a tree via the
    `Read` subobject attached this instance."""

class Tree(handlers.Handler):
    """Tree Handler. Data is m

    Parameters
    ----------
    tree : string
        a tree in the form of a newick string.
    """
    def __init__(self, tree, **kwargs):
        super(Tree, self).__init__(tree=tree, **kwarg)
        self._dendropy = Dendropy(schema="")

    @property
    def dendropy(self):
        try:
            return self._dendropy
        except AttributeError:
            raise DendropyError

    def link_homologs(self):
        """Link homolog handlers to tree.
        """

    def link_ancestors(self):
        """Link ancestor handlers to tree.
        """

    def remove(self):
        """Remove Homologs from tree.
        """

    def reroot(self):
        """Reroot the tree on a given node or branch.
        """

    def subtree(self):
        """Return a subtree from the main tree.
        """

    def prune(self):
        """Remove a homolog, ancestor (and children), or branch from tree.
        """

    def place(self):
        """Add new subtree into tree.
        """

class TreeList(handlers.HandlerContainer):
    """Container for Trees.
    """
    def __init__(self, *Trees, **kwargs):
        super(TreeList, self).__init__(*Trees, **kwargs)

    @property
    def _prefix(self):
        return "Tree"

    @property
    def _child_type(self):
        """The Handler object subclass passed into this Container. Must be
        a contained in a list. Can take in many types of Handlers
        """
        return [Tree]
