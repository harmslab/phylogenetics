from . import handlers

class Tree(handlers.Handler):
    """
    """
    def __init__(self, **kwargs):
        super(Tree, self).__init__(**kwarg)

    def link_homologs(self):
        """
        """

    def link_ancestors(self):
        """
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

class TreeList(handlers.HandlerContainer)
