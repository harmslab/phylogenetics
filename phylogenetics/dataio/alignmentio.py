# Module for input and output of alignment.

from . import base
from .formats import fasta, phylip, json, pickle

class Write(base.Write):

    def __init__(self, Alignment):
        """ Object for reading and writing alignment data from a homolog set. """
        self._Alignment = Alignment

    def _object_to_sequences(self, alignment="latest", tags=["id"]):
        """ Convert alignment data to sequence data type. """
        sequences = []
        alignment_dict = self._Alignment._alignments[alignment]
        homologs = self._Alignment._HomologSet.homologs.items()
        for id, homolog in homologs:
            data = homolog.get(*tags)
            sequences.append((tuple([data[t] for t in tags]), alignment_dict[id]))
        return sequences

    def _object_to_data(self):
        """ Write alignment data to metadata format.

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
        return self._Alignment._alignments

    @base.write_to_file
    def fasta(self, alignment="latest", tags=["id"]):
        """ Write alignment to fasta format. """
        data = self._object_to_sequences(alignment=alignment, tags=tags)
        output = fasta.write(data)
        return output

    @base.write_to_file
    def pickle(self, alignment="latest"):
        """Write alignment to pickle string.
        """
        data = self._object_to_data(alignment=alignment)
        output = pickle.write(data)
        return output

    @base.write_to_file
    def json(self, alignment="latest"):
        """Write alignment to json string.
        """
        data = self._object_to_data(alignment=alignment)
        output = json.write(data)
        return output

    @base.write_to_file
    def phylip(self, alignment="latest"):
        """" Write alignment to phylip format. """
        data = self._object_to_sequences(alignment=alignment)
        output = phylip.write(data)
        return output

    @base.write_to_file
    def csv(self, tags=[]):
        """ Write alignment to csv format."""
        sequence_metadata = self._object_to_data(tags=None)
        output = csv.write(sequence_metadata)
        return output

class Read(base.Read):

    def __init__(self, Alignment):
        """ Module for reading alignment for Homologset from various formats. """
        self._Alignment = Alignment

    def _sequences_to_object(self, data, tags=["id"]):
        """ Add sequence_data to alignment

        example:
        data = [
            (("XX00000001", "dog"), "ASHASHSAEFASHAS"),
            (("XX00000002", "cat"), "ASTASHSAASDGAWE"),
            ...
        ]


        """
        # Check for ids in alignment
        try:
            index = tags.index("id")
        except ValueError:
            raise Exception(""" One tag in alignment must be `id`. """)

        # Move old alignment to prevent overwriting.
        self._Alignment._move_old_alignment()

        # Initialize a list to old all ids in file.
        ids_in_alignment_file = list()

        alignment = {}
        # Iterate through the sequence data
        for pair in data:
            attributes = pair[0]
            sequence = pair[1]

            id = attributes[index]

            if len(attributes) != len(tags):
                raise Exception(""" Number of attributes do not match number of tags given. """)

            alignment[id] = sequence
            # Track ids that are in alignment
            ids_in_alignment_file.append(id)

        # Set this alignment as the latest.
        self._Alignment._alignments["latest"] = alignment

        # Check if Homologs were removed from Alignment file. If so, remove from
        # HomologSet.
        ids_in_alignment_file = set(ids_in_alignment_file)
        ids_in_HomologSet = set(self._Alignment._HomologSet.list_ids)
        diff = list(ids_in_HomologSet - ids_in_alignment_file)

        self._Alignment._HomologSet.rm(diff)

        return self._Alignment

    def _data_to_object(self, data):
        """Add phylip data to alignment.

        Parameters
        ----------

        Returns
        -------

        """
        # Initialize a list to old all ids in file.
        ids_in_alignment_file = [homolog["id"] for homolog in data]
        # Move old alignment to prevent overwriting.
        self._Alignment._move_old_alignment()
        # Set this alignment as the latest.
        self._Alignment._alignments["latest"] = data
        # Check if Homologs were removed from Alignment file. If so, remove from
        # HomologSet.
        ids_in_alignment_file = set(ids_in_alignment_file)
        ids_in_HomologSet = set(self._Alignment._HomologSet.list_ids)
        diff = list(ids_in_HomologSet - ids_in_alignment_file)
        self._Alignment._HomologSet.rm(diff)
        return self._Alignment

    @base.read_from_file
    def fasta(self, data, tags=["id"]):
        """ Read alignment from fasta file. """
        sequence_data = fasta.read(data, tags=tags)
        self._data_to_object(sequence_data)
        return self._Alignment

    @base.read_from_file
    def phylip(self, data):
        """ Read alignment from fasta file. """
        sequence_data = phylip.read(data)
        self._sequences_to_object(sequence_data)
        return self._Alignment

    @base.read_from_file
    def json(self, data):
        """Read json data and add to alignment.
        """
        metadata = json.read(data)
        self._data_to_object(metadata)
        return self._Alignment

    @base.read_from_file
    def pickle(self, data):
        """Read pickle data and add to alignment.
        """
        metadata = pickle.read(data)
        self._data_to_object(metadata)
        return self._Alignment
