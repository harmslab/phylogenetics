from .base import read_from_file, write_to_file


# import objects to bind to Project class
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
                        rst)

class Read(object):

    def __init__(self, Project):
        """"""
        self._Project = Project

    @read_from_file
    def fasta(self, data, tags=("id",)):
        """Read fasta file.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.fasta(data, tags=tags)

    @read_from_file
    def alignment(self, data, tags=("id",)):
        """Read alignment from file.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.add( Alignment(self._Project.HomologSet) )
        self._Project.Alignment.Read.fasta(data, tags=tags)

    @read_from_file
    def entre_xml(self, data):
        """Read downloaded XML from Entrez server.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.HomologSet.Read.entre_xml(data)

    @read_from_file
    def phylip(self, data):
        """Read phylip file.
        """
        if hasattr(self._Project, "HomologSet") is False:
            # Add a HomologSet to project
            self._Project.add( HomologSet() )

        self._Project.add( Alignment(self._Project.HomologSet) )
        self._Project.Alignment.Read.phylip(data, tags=tags)

    @read_from_file
    def csv(self, data):
        """Read csv file.
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
        """
        
