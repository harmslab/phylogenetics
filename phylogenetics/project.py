import importlib

from . import handlers
from . import sequences
from . import alignments
from . import trees
from . import ancestors

from .dataio import read, write, formats
from .exttools import msaprobs, phyml, paml

def checkpoint(method):
    """Wrap an project's method tp save a checkpoint of the object after performing
    the function
    """
    def write_checkpoint(*args, **kwargs):
        """
        """
        project = args[0]
        return method(*args, **kwargs)
    return checkpoint

class BadHandlerError(Exception):
    """Error message for passing in a bad handler to the Project class."""

class Project(handlers.HandlerContainer):
    """Main object for managing a phylogenetics project.
    """
    def __init__(self, **kwargs):
        super(Project, self).__init__(**kwargs)

    @property
    def _schemas(self):
        """Read/Write schemas available to project class.
        """
        return ["pickle", "json"]

    @property
    def _child_types(self):
        return [sequences.SequenceList,
            alignments.AlignmentList,
            trees.TreeList,
            ancestors.AncestorList
        ]

    def add(self, *Handlers):
        """Add phylogenetic handlers to a Project class.
        """
        # Add items to the project class
        for Handler in Handlers:
            itemname = Handler.type
            # Try to add the Handler. If its a bad Handler, raise exception.
            try:
                method = getattr(self, "_add_" + itemname)
                method(Handler)
            except:
                raise BadHandlerError(itemname + " is not a valid object for Project class.")

    def link(self, **links):
        """Create a link between two objects.
        """
        # Add manual links between objects
        for source, target in links.items():
            Source = getattr(self, source)
            Target = getattr(self, target)
            Source.link(Target)

    def unlink(self):
        """"""

    @handlers.history
    def _align(self, **options):
        """
        """
        self.SequenceList
        output = msaprobs.run(**options)
        alignment = alignments.Alignment()
        alignment.read(path=output, schema="fasta")


    @handlers.history
    def _tree(self):
        """"""

    @handlers.history
    def _reconstruct(self):
        """"""

    def _add_SequenceList(self, SequenceList):
        """Add HomologSet Set to PhylogeneticsProject object."""
        # Expose the align method of this object to user
        self.SequenceList=SequenceList
        self._contents["SequenceList"] = self.SequenceList
        setattr(self, "align", self._align)

    def _add_AlignmentList(self, AlignmentList=None):
        """Initialize an AlignmentList for the project.
        """
        if AlignmentList is None:
            AlignmentList=alignments.AlignmentList()
        self.AlignmentList=AlignmentList
        self._contents["AlignmentList"] = self.AlignmentList

    def _add_Alignment(self, Alignment):
        """Add Alignment to PhylogeneticsProject object."""
        try:
            self.AlignmentList.add(Alignment)
        except AttributeError:
            self._add_AlignmentList()
            self.AlignmentList.add(Alignment)
        setattr(self, "tree", self._tree)

    def _add_TreeList(self, TreeList=None):
        """Initialize an TreeList for the project.
        """
        if TreeList is None:
            TreeList=trees.TreeList()
        self.TreeList=TreeList
        self._contents["TreeList"] = self.TreeList

    def _add_Tree(self, Tree):
        """Add Tree to Phylogenetics Project object."""
        try:
            self.TreeList.add(Tree)
        except AttributeError:
            self._add_TreeList()
            self.TreeList.add(Tree)
        # Expose the reconstruction methods of this project object
        setattr(self, "reconstruct", self._reconstruct)

    def _add_AncestorList(self):
        """Add a AncestorSet object to PhylogeneticsProject object."""
        self.AncestorList=AncestorList
        self._contents["AncestorList"] = self.AncestorList
