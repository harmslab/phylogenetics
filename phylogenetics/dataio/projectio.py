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

    def _data_to_object(data):
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

    def _data_to_sequences(data):
        """This method doesn't do anything.
        """
        print("""Use internal objects (i.e. HomologSet, Alignment) to read \
        sequences from file, or the `files method`.""")

    @base.read_from_file
    def pickle(self, data):
        """Read a pickled string for project metadata."""
        metadata = pickle.read(data)
        self._data_to_object(metadata)
        return self._Project

    @base.read_from_file
    def json(self, data):
        """Read a json string for project metadata."""
        metadata = json.read(data)
        self._data_to_object(metadata)
        return self._Project

class Write(base.Write):
    """Writing object for project class.
    """
    def __init__(self, Project):
        self._Project = Project

    def _object_to_data(self):
        """Write project object to metadata dictionary."""
        data = {}
        for key, object self._Project._components.items():
            data[key] = object._object_to_data()
        return data

    def _object_to_sequences(self):
        """This method doesn't do anything in this object."""
        print("""Use internal objects (i.e. HomologSet, Alignment) to write \
        sequence data to file.""")

    @base.write_to_file
    def pickle(self):
        """Write phylogenetics project to pickle string."""
        metadata = self._object_to_data()
        pickle.write(metadata)

    @base.write_to_file
    def json(self):
        """Write phylogenetics project to json string."""
        metadata = self._object_to_data()
        json.write(metadata)
