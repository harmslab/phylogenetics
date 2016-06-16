"""Module for reading input and writing output of Homolog Object.

Writes homolog object attributes out as a python dictionary with metadata which can
easily be written:


Example :
>>> { "XX00000000" :
        {
            "sequence" : "XXXXXX ....",
            "latest_align" : "XX----XX-X--XX...",
            "organism" : "mouse"
        }
    }
"""

# ----------------------------------------------
# Imports
# ----------------------------------------------

from . import base
from .formats import fasta, csv, json, pickle

# ----------------------------------------------
# Winting moduel
# ----------------------------------------------

class Write(base.Write):
    """ Writing object for Homolog Object.

    """
    def __init__(self, Homolog):
        self._Homolog = Homolog

    def _object_to_sequences(self, tags=[]):
        """ Write Homolog as sequence_data datatype.

        Arguments
        ---------
        tags : list
            list of attributes in Homolog object to write

        Returns
        -------
        sequence_data : dict

        """
        # Built a tuple of tags
        tag_data = tuple()
        for t in tags:
            tag_data += (getattr(self._Homolog, t),)

        # Built tuple pair of tags to sequence
        sequence_data = (tag_data, self._Homolog.sequence)

        return sequence_data

    def _object_to_data(self):
        """
        """
        return self._Homolog.attrs

    @base.write_to_file
    def fasta(self, tags=[], aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        sequence_data = self._object_to_sequences()
        output = fasta.write(sequence_data)
        return output

    @base.write_to_file
    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.write(self._object_to_data())

    @base.write_to_file
    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.write(self._object_to_data())

    @base.write_to_file
    def csv(self, tags=None, header=True, delimiter=",", **kwargs):
        """ write csv string. """
        # Get all attributes if tags are not specified
        if tags is None:
            tags = list(vars(self._Homolog).keys())

        # Add header is specified
        if header:
            f = delimiter.join(tags)
        else:
            f = ""

        # Try to print tag. If it doesnt exist, then leave blank
        vals = list()
        for t in tags:
            try:
                vals.append(str(getattr(self._Homolog,t)))
            except:
                vals.append("")
        f += delimiter.join(vals) + "\n"

        return f

class Read(base.Read):

    def __init__(self, Homolog):
        """ Object for reading and writing HomologSets. """
        self._Homolog = Homolog

    def _sequences_to_object(self, sequence_data, tags=None):
        """ Method to take transform sequence data-structure from reading methods
            to a Homolog.
        """
        # Split the sequence tuple into relevent parts
        attributes = sequence_data[0]
        sequence = sequence_data[1]

        # Fill in tags if none are given
        if tags is None:
            # Either use unique ID or set as defline
            if attributes[0][0:2] == "XX" and len(attributes[0]) == 10:
                tag = "id"
            else:
                tag = "defline"

            self._Homolog.addattr(tag, "|".join(attributes))
        else:
            # Add attributes in the header
            for i in len(tags):
                self._Homolog.addattr(tags[i], attributes[i])

        # Add sequence to Homolog object
        self._Homolog.addattr("sequence", sequence)

        return self._Homolog

    def _data_to_object(self, sequence_metadata):
        """Use read take to populate a Homolog object.

        Arguments
        ---------
        data : dict
            dict with all the metdata.
        """
        for key, value in sequence_metadata.items():
            self._Homolog.addattr(key, value)
        return self._Homolog

    @base.read_from_file
    def fasta(self, data, tags=None):
        """ Read a fasta string.

            Arguments:
            ---------
            data: None

            tag: tuple
                attributes to add fasta headers.

        """
        # Read fasta file
        sequence_data = fasta.read(data)
        return self._sequences_to_object(sequence_data)

    @base.read_from_file
    def csv(self, data):
        pass

    @base.read_from_file
    def json(self, data):
        pass
