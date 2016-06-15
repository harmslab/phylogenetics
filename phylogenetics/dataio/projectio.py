"""

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

    ""



}
"""
from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.alignment import Alignment
from phylogenetics.tree import Tree
from phylogenetics.ancestors import Ancestor, AncestorSet
from phylogenetics.reconstruction import Reconstruction
from .base import read_from_file, write_to_file

from .formats import (csv,
                        entrez_xml,
                        fasta,
                        phylip,
                        pickle,
                        rst)

class Read(object):
    """Module for reading in files to project object.
    """
    def __init__(self, Project):
        self._Project = Project

    @read_from_file
    def fasta(self, data, tags=("id",)):
        """Read fasta string.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.fasta(data, tags=tags)

    @read_from_file
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

    @read_from_file
    def entrez_xml(self, data):
        """Read downloaded XML string from Entrez server.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.entre_xml(data)

    @read_from_file
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

    @read_from_file
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

    @read_from_file
    def rst(self, data):
        """Read PAML output form file.

        Note: to read from a file, relace data argument with `fname=` keyword
        argument.
        """
