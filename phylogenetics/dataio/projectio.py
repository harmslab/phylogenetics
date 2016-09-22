from __future__ import absolute_import
from collections import OrderedDict

from . import base
from phylogenetics.sequences import Sequence, SequenceLists
from phylogenetics.alignment import Alignment, AlignmentList
from phylogenetics.tree import Tree, TreeList
from phylogenetics.ancestors import Ancestor, AncestorSet, AncestorSetList

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

    def _data_to_object(self, data):
        """Reads a project object from dictionary metadata.
        """
        objects = {
            #"Alignment" : Alignment,
            "AlignmentList" : AlignmentList,
            #"Ancestor" : Ancestor,
            #"AncestorSet" : AncestorSet,
            "AncestorSetList" : AncestorSetList,
            #"Sequence" : Sequence
            "SequenceList" : SequenceList,
            #"Tree" : Tree,
            "TreeList" : TreeList
        }

        for 

        for name, item in meta.items():
            if name in data:
                key = name
                # Initialize the object
                new_object = meta[key]
                # Read the data
                new_object.Read._data_to_object(data[key])
                # Add to project object
                self._Project.add(new_object)

    def _data_to_sequences(self, data):
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
        for key, obj in self._Project._components.items():
            data[key] = obj.Write._object_to_data()
        return data

    def _object_to_sequences(self):
        """This method doesn't do anything in this object."""
        print("""Use internal objects (i.e. HomologSet, Alignment) to write \
        sequence data to file.""")

    @base.write_to_file
    def pickle(self):
        """Write phylogenetics project to pickle string."""
        data = self._object_to_data()
        metadata = pickle.write(data)
        return metadata

    @base.write_to_file
    def json(self):
        """Write phylogenetics project to json string."""
        data = self._object_to_data()
        metadata = json.write(data)
        return metadata
