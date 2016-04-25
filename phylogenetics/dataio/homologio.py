# Module for input/output of Homolog Object

from .base import read_from_file, write_to_file
from .formats import fasta, csv

class Write(object):

    def __init__(self, Homolog):
        """ Writing object for Homolog Object. """
        self._Homolog = Homolog

    def _homolog_to_sequence_data(self, tags=("id",)):
        """ Write Homolog as sequence_data datatype.

            Arguments:
            ---------
            sequence_data


            Output Format:
            -------------

        """
        # Built a tuple of tags
        tag_data = tuple()
        for t in tags:
            tag_data += (getattr(self._Homolog, t),)

        # Built tuple pair of tags to sequence
        sequence_data = (tag_data, self._Homolog.sequence)

        return sequence_data

    def _homolog_to_sequence_metadata(self, tags=None):
        """ Write Homolog to sequence_metadata datatype. """
        sequence_metadata = []
        if tags is None:
            sequence_metadata.append(self._Homolog.attrs)
        else:
            metadata = dict([(t, getattr(self._Homolog, t)) for t in tags])
            sequence_metadata.append(metadata)
        return sequence_metadata


    @write_to_file
    def fasta(self, tags=("id",), aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        sequence_data = self._homolog_to_sequence_data()
        output = fasta.write(sequence_data)
        return output

    @write_to_file
    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.dumps(self._Homolog.attrs)


    @write_to_file
    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self._Homolog)

    @write_to_file
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

class Read(object):

    def __init__(self, Homolog):
        """ Object for reading and writing HomologSets. """
        self._Homolog = Homolog

    def _sequence_data_to_homolog(self, sequence_data, tags=None):
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

    def _sequence_metadata_to_homolog(self, sequence_metadata):
        """ Convert a sequence metadata datatype to a homolog object. """
        for key, value in sequence_metadata.items():
            self._Homolog.addattr(key, value)
        return self._Homolog

    @read_from_file
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
        return self._sequence_data_to_homolog(sequence_data)

    @read_from_file
    def csv(self, data):
        pass

    @read_from_file
    def json(self, data):
        pass
