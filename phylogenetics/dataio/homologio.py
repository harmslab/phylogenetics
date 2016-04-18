# Module for input/output of Homolog Object

from .base import read_from_file, write_to_file
from phylogenetics.dataio import fasta

class Write(object):

    def __init__(self, Homolog):
        """ Writing object for Homolog Object. """
        self._Homolog = Homolog

    def _homolog_to_sequence_data(self, tags=("id",)):
        """ Write homolog set as sequence data to be written to file.

            Arguments:
            ---------
            sequence_data
        """
        # Built a tuple of tags
        tag_data = tuple()
        for t in tags:
            tag_data += (self.Homolog.get_attr(t),)

        # Built tuple pair of tags to sequence
        sequence_data = (tag_data, Homolog.sequence)

        return sequence_data

    @write_to_file
    def fasta(self, tags=("id",), aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        sequence_data = self._homolog_to_sequence_data()
        fasta.write()
        return f

    @write_to_file
    def phylip(self, alignment_index=None, **kwargs):
        """ Return a PhyML formatted string to write to file.

            Only allowed when Homolog has alignment attribute
        """
        # Check that an alignment attribute exists in homolog.
        if hasattr(self._homolog, "latest_align") is False:
            raise Exception("Homolog must have an alignment attribute to write phylip.")

        # Use specified alignment key
        if tags is not None:
            alignment = getattr(self._homolog, tags)
        else:
            alignment = self._homolog.latest_align

        f = "%s\n%s\n" % (self._homolog.id, alignment)
        return f


    @write_to_file
    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.dumps(self._homolog.attrs)


    @write_to_file
    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self._homolog)


    @write_to_file
    def csv(self, tags=None, header=True, delimiter=",", **kwargs):
        """ write csv string. """
        # Get all attributes if tags are not specified
        if tags is None:
            tags = list(vars(self._homolog).keys())

        # Add header is specified
        if header:
            f = delimiter.join(tags)
        else:
            f = ""

        # Try to print tag. If it doesnt exist, then leave blank
        vals = list()
        for t in tags:
            try:
                vals.append(str(getattr(self._homolog,t)))
            except:
                vals.append("")
        f += delimiter.join(vals) + "\n"

        return f

class Read(object):

    def __init__(self, Homolog):
        """ Object for reading and writing HomologSets. """
        self._Homolog = Homolog

    def _sequence_data_to_homolog(self, sequence_data):
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

    @read_from_file
    def phylip(self, data):
        pass
