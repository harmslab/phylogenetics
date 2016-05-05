# Module for input and output of HomologSet objects
from phylogenetics import homologs
from .base import write_to_file, read_from_file
from .formats import fasta, pickle, csv, entrez_xml

class Write(object):
    """Object for writing out the metadata of a HomologSet object.

    There are two basic output data-structures for the HomologSet.

    1. Sequence data::

        data = [
            (("XX00000001", "dog"), "ASHASHSAEFASHAS"),
            (("XX00000002", "cat"), "ASTASHSAASDGAWE"),
            ...
        ]

    2. Sequence Metadata::

        sequence_metadata = [
            {
                "species": "seq0",
                "organism": "agasda",
                "sequence": "SHDAHADJAEAHASDASDHASDGBSHERW",
            },
            {
                "species": "seq1",
                "organism": "afher",
                "sequence": "OENGBSDMLWETJALSGMSDALGMASDFW",
            },
            ...
        ]

    """
    def __init__(self, HomologSet):
        self._HomologSet = HomologSet

    def _homologset_to_sequence_data(self, tags=("id",)):
        """ Convert alignment data to sequence data type """
        sequence_data = list()
        for id, homolog in self._HomologSet.homologs.items():
            # Construct tag list
            labels = tuple()

            for t in tags:
                labels += (getattr(homolog, t),)
            # Append tuple pair to list
            sequence_data.append((labels, getattr(homolog, "sequence")))

        return sequence_data

    def _homologset_to_metadata(self, tags=("id",)):
        """ Write alignment data to metadata format.
        """
        sequence_metadata = []
        if tags is None:
            for id, homolog in self._HomologSet.homologs.items():
                sequence_metadata.append(homolog.attrs)
        else:
            for id, homolog in self._HomologSet.homologs.items():
                metadata = dict([(t, getattr(homolog, t)) for t in tags])
                sequence_metadata.append(metadata)
        return sequence_metadata

    @write_to_file
    def fasta(self, tags=("id",)):
        """ Return string in fasta format for the set."""
        sequence_data = self._homologset_to_sequence_data(tags=tags)
        output = fasta.write(sequence_data)
        return output

    @write_to_file
    def pickle(self):
        """ Return pickle string of homolog set. """
        return pickle.write(self._HomologSet)

    @write_to_file
    def csv(self, tags=("id",), delimiter=","):
        """ Return csv string. """
        sequence_metadata = self._homologset_to_metadata(tags=tags)
        output = csv.write(sequence_metadata, delimiter=delimiter)
        return output


class Read(object):
    """Object bound to HomologSets for reading homolog set data.

    There are two basic output data-structures for the HomologSet.

    1. Sequence data::

        data = [
            (("XX00000001", "dog"), "ASHASHSAEFASHAS"),
            (("XX00000002", "cat"), "ASTASHSAASDGAWE"),
            ...
        ]

    2. Sequence Metadata::

        sequence_metadata = [
            {
                "species": "seq0",
                "organism": "agasda",
                "sequence": "SHDAHADJAEAHASDASDHASDGBSHERW",
            },
            {
                "species": "seq1",
                "organism": "afher",
                "sequence": "OENGBSDMLWETJALSGMSDALGMASDFW",
            },
            ...
        ]

    Parameters
    ----------
    HomologSet : phylogenetics.homolog.HomologSet object
        HomologSet to which the reading object will be bound.

    Examples
    --------

    """
    def __init__(self, HomologSet):
        self._HomologSet = HomologSet

    def _sequence_data_to_homologset(self, data, tags=("id",)):
        """ Add sequence_data to alignment
        """
        # Check for ids in alignment
        try:
            index = tags.index("id")
        except ValueError:
            raise Exception(""" One tag in alignment must be `id`. """)

        # Iterate through the sequence data
        for pair in data:
            attributes = pair[0]
            sequence = pair[1]

            # Get ID from attributes list.
            id = attributes[index]

            # Check that attribute list and tag list are the same length.
            if len(attributes) != len(tags):
                raise Exception(""" Number of attributes do not match number of tags given. """)

            # Update homologs that are already in the Set.
            if hasattr(self._HomologSet, id):
                # Get Homolog object
                Homolog = getattr(self._HomologSet, id)
            # otherwise, add a new Homolog object
            else:
                Homolog = homologs.Homolog(id)
                self._HomologSet.add(Homolog)

            # Add updated attributes
            for i in range(len(tags)):
                Homolog.addattr(tags[i], attributes[i])

            Homolog.addattr("sequence", sequence)

        return self._HomologSet


    def _sequence_metadata_to_homologset(self, sequence_metadata):
        """ Add sequence_metadata to alignment.
        """
        for s in sequence_metadata:
            # Update homologs that are already in the Set.
            if "id" in s:
                # Get Homolog object
                Homolog = getattr(self._HomologSet, s["id"])
            # otherwise, add a new Homolog object
            else:
                id = "XX%08d" % int(self._HomologSet.max_id + 1)
                Homolog = homologs.Homolog(id)
                self._HomologSet.add(Homolog)

            # Iterate over attributes in s and add to homolog.
            for key, value in s.items():
                Homolog.addattr(key, value)

        return self._HomologSet

    @read_from_file
    def fasta(self, data, tags=("id",)):
        """ Add sequence data from fasta to HomologSet. """
        sequence_data = fasta.read(data)
        self._sequence_data_to_homologset(sequence_data, tags=tags)
        return self._HomologSet

    @read_from_file
    def csv(self, data):
        """ Add csv data to HomologSet. """
        sequence_metadata = csv.read(data)
        self._sequence_metadata_to_homologset(sequence_metadata)
        return self._HomologSet

    @read_from_file
    def entrez_xml(self, data):
        """ Add entrez xml output to homolog. """
        sequence_metadata = entrez_xml.read(data)
        self._sequence_metadata_to_homologset(sequence_metadata)
        return self._HomologSet
