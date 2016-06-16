# Module for input/output of Homolog Object

import .base
from .formats import fasta, csv, rst

class Write(base.Write):

    def __init__(self, Ancestor):
        """ Writing object for Homolog Object. """
        self._Ancestor = Ancestor

    def _ancestor_to_sequence_data(self, tags=("id",)):
        """ Write Ancestor as sequence_data datatype.

            Arguments:
            ---------
            sequence_data


            Output Format:
            -------------

        """
        # Built a tuple of tags
        tag_data = tuple()
        for t in tags:
            tag_data += (getattr(self._Ancestor, t),)

        # Built tuple pair of tags to sequence
        sequence_data = (tag_data, self._Ancestor.sequence)

        return sequence_data

    def _ancestor_to_sequence_metadata(self, tags=None):
        """ Write Ancestor to sequence_metadata datatype. """
        sequence_metadata = []
        if tags is None:
            sequence_metadata.append(self._Ancestor.attrs)
        else:
            metadata = dict([(t, getattr(self._Ancestor, t)) for t in tags])
            sequence_metadata.append(metadata)
        return sequence_metadata


    @base.write_to_file
    def fasta(self, tags=("id",), aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        sequence_data = self._Ancestor_to_sequence_data()
        output = fasta.write(sequence_data)
        return output

    @base.write_to_file
    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.dumps(self._Ancestor.attrs)


    @base.write_to_file
    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self._Ancestor)

    @base.write_to_file
    def csv(self, tags=None, header=True, delimiter=",", **kwargs):
        """ write csv string. """
        # Get all attributes if tags are not specified
        if tags is None:
            tags = list(vars(self._Ancestor).keys())

        # Add header is specified
        if header:
            f = delimiter.join(tags)
        else:
            f = ""

        # Try to print tag. If it doesnt exist, then leave blank
        vals = list()
        for t in tags:
            try:
                vals.append(str(getattr(self._Ancestor,t)))
            except:
                vals.append("")
        f += delimiter.join(vals) + "\n"

        return f

class Read(base.Read):

    def __init__(self, Ancestor):
        """ Object for reading and writing Ancestor. """
        self._Ancestor = Ancestor

    def _ancestor_data_to_ancestor(self, data):
        """ Method to take transform sequence data-structure from reading methods
            to a Ancestor.
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

            self._Ancestor.addattr(tag, "|".join(attributes))
        else:
            # Add attributes in the header
            for i in len(tags):
                self._Ancestor.addattr(tags[i], attributes[i])

        # Add sequence to Ancestor object
        self._Ancestor.addattr("sequence", sequence)

        return self._Ancestor

    def _sequence_metadata_to_ancestor(self, sequence_metadata):
        """ Convert a sequence metadata datatype to a ancestor object. """
        for key, value in sequence_metadata.items():
            self._Ancestor.addattr(key, value)
        return self._Ancestor

    @base.read_from_file
    def rst (self, data):
        """ Read PAMLS's output. """
        pass

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
        return self._sequence_data_to_ancestor(sequence_data)

    @base.read_from_file
    def csv(self, data):
        pass

    @base.read_from_file
    def json(self, data):
        pass
