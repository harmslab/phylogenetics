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


def project_metadata(schema, metadata):
    """Implant sub-metadata into a project metadata dictionary.
    """
    options = {
        "fasta" : dict(
            type="Project",
            module="phylogenetics.project",
            contents=[metadata]
        ),
        "ali" : dict(
            type="Project",
            module="phylogenetics.project",
            contents=[dict(
                type="AlignmentList",
                module="phylogenetics.alignments",
                contents=[metadata]
            )]
        ),
        "newick" : dict(
            type="Project",
            module="phylogenetics.project",
            contents=[dict(
                type="TreeList",
                module="phylogenetics.trees",
                contents=[metadata]
            )]
        ),
        "rst" : dict(
            type="Project",
            module="phylogenetics.project",
            contents=[dict(
                type="AncestorList",
                module="phylogenetics.ancestors",
                contents=[metadata]
            )]
        ),
    }
    return options[schema]


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
        return ["pickle", "json", "fasta", "ali"]

    @property
    def _child_types(self):
        return [sequences.SequenceList,
            alignments.AlignmentList,
            trees.TreeList,
            ancestors.AncestorList
        ]

    @handlers.history
    @read.file
    def read(self, data, schema):
        """Read data with schema to object."""
        parser = getattr(formats, schema).read
        metadata = parser(data)
        metadata = project_metadata(schema, metadata)
        self._metadata_to_object(metadata)

    def _add(self, *Handlers):
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
    def _align(self, sequence_list_name="SequenceList", **options):
        """Run alignment using MSAProbs.
        """
        fname = "pre-alignment"
        path = fname + ".fasta"
        # Get list of sequences to align
        SequenceList = getattr(self, sequence_list_name)
        # Write sequence list to disk
        SequenceList.write(path=path, schema="fasta")
        # Run alignment and read it
        output = msaprobs.run(fasta_fname=fname,**options)
        alignment = alignments.Alignment()
        alignment.read(path=output, schema="fasta")
        # Add to project
        self._add_Alignment(alignment)

    @handlers.history
    def _tree(self):
        """"""

    @handlers.history
    def _reconstruct(self):
        """"""

    def _add_SequenceList(self, SequenceList):
        """Add HomologSet Set to PhylogeneticsProject object."""
        # Expose the align method of this object to user
        # Enforce that the id is standardized and not numbered
        SequenceList._addattr(id="SequenceList")
        setattr(self, SequenceList.id, SequenceList)
        self._contents[SequenceList.id] = SequenceList
        setattr(self, "align", self._align)

    def _add_AlignmentList(self, AlignmentList=None):
        """Initialize an AlignmentList for the project.
        """
        AlignmentList._addattr(id="AlignmentList")
        if AlignmentList is None:
            AlignmentList=alignments.AlignmentList()
        setattr(self, AlignmentList.id, AlignmentList)
        self._contents[AlignmentList.id] = AlignmentList

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
        TreeList._add(id="TreeList")
        if TreeList is None:
            TreeList=trees.TreeList()
        setattr(self, TreeList.id, TreeList)
        self._contents[TreeList.id] = TreeList

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
