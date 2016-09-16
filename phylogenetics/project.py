from . import handlers
from . import sequences
from . import alignments
from . import trees
from . import ancestors

def create_checkpoint(method):
    """Wrap an object's method save a checkpoint of the object before performing
    the function
    """
    def checkpoint(*args, **kwargs):
        return method(*args, **kwargs)
    return checkpoint

class Project(handlers.HandlerContainer):
    """Main object for managing a phylogenetics project.
    """
    def __init__(self, *args, **kwargs):
        self._child_types_ = [sequences.SequenceList,
            alignments.AlignmentList,
            trees.TreeList,
            ancestors.AncestorList
        ]
        super(Project, self).__init__(*args, **kwargs)

    @property
    def _prefix(self):
        return "Project"

    @property
    def _child_types(self):
        """The Handler object subclass passed into this Container. Must be
        a contained in a list. Can take in many types of Handlers
        """
        return self._child_types_

    def add_child_type(self, obj):
        """Add a child type to project"""
        self._child_type_.append(obj)

    def save(self, fname):
        """Write metadata of handler container to json file.
        """
        with open(fname, "w") as f:
            json.dump(self.metdata, f)

    def load(self, fname):
        """Load in metadata from a project.
        """
        with open(fname, "r") as f:
            data = json.load(f)

#        for key, value in data.items():
