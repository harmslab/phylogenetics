"""Input/Output module for saving phylogenetics Project objects to disk.

Writes metadata as so:

{
    "HomologSet" : [
        {
            "id" : "XX00000000",
            "sequence" : "XJSMHTELTWF...",
            "accver" : "#########",
        },
        {
            "id" : "XX00000001",
            "sequence" : "DSHASTEJOGASHJ...",
            "accver" : "######"
        }
    ],
    "Alignment" :

    "Tree" : (,),

    "AncestorSet" : [
        {
            "id" : "ZZ00000000",
            "sequence" : "XJSMHTELTWF...",
            "posterior" : "0.99"
        },
        {
            "id" : "ZZ00000001",
            "sequence" : "DSHASTEJOGASHJ...",
            "posterior" : "0.95"
        }
    ]
}
"""

import .base
from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.alignment import Alignment
from phylogenetics.tree import Tree
from phylogenetics.ancestors import Ancestor, AncestorSet
from phylogenetics.reconstruction import Reconstruction

from .formats import (csv,
                        entrez_xml,
                        fasta,
                        phylip,
                        pickle,
                        json,
                        rst)

class Read(base.Read):
    """Module for reading in files to project object.
    """
    def __init__(self, Project):
        self._Project = Project

    def _data_to_project(data):
        """Reads a project object from dictionary metadata.
        """
        object_types = {
            "HomologSet" : HomologSet,
            "Alignment" : Alignment,
            "Tree" : Tree,
            "AncestorSet" : AncestorSet,
            "Reconstruction" : Reconstruction,
        }

        # walk through data object and add each feature to project class.
        for key in data:
            new_object = object_types[key]
            new_object.Read._data_to_object(data[key])
            self._Project.add(new_object)

    def _sequences_to_project(data):
        """Read project object from sequence data.
        """
        pass


    @base.read_from_file
    def fasta(self, data, tags=("id",)):
        """Read fasta string.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.fasta(data, tags=tags)

    @base.read_from_file
    def alignment(self, data, tags=("id",)):
        """Read alignment from string.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.add( Alignment(self._Project.HomologSet) )
        self._Project.Alignment.Read.fasta(data, tags=tags)

    @base.read_from_file
    def entrez_xml(self, data):
        """Read downloaded XML string from Entrez server.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.entre_xml(data)

    @base.read_from_file
    def phylip(self, data):
        """Read phylip string.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.add( Alignment(self._Project.HomologSet) )
        self._Project.Alignment.Read.phylip(data, tags=tags)

    @base.read_from_file
    def csv(self, data):
        """Read csv string.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.csv(data)

        # Check if an alignment is in the csv.
        if hasattr(list(project.HomologSet.homologs.values())[0], "latest_align"):
            # Add alignment to project
            self._Project.add( Alignment( self._Project.HomologSet ) )

    @base.read_from_file
    def rst(self, data):
        """Read PAML output form file.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """

class Write(base.Write):
    """Writing object for project class.
    """
    def __init__(self, Project):
        self._Project = Project
