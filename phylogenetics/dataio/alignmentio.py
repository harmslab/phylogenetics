# Module for input and output of alignment.

from .base import read_from_file, write_to_file
from .formats import fasta, phylip

class Write(object):

    def __init__(self, Alignment):
        """ Object for reading and writing alignment data from a homolog set. """
        self._Alignment = Alignment

    def _align_to_sequence_data(self, alignment="latest_align", tags=("id",)):
        """ Convert alignment data to sequence data type """
        sequence_data = list()
        for id, homolog in self._Alignment._HomologSet.homologs.items():
            # Construct tag list
            labels = tuple()
            for t in tags:
                labels += (getattr(homolog, t),)
            # Append tuple pair to list
            sequence_data.append((labels, getattr(homolog, alignment)))

        return sequence_data

    def _align_to_phylip_data(self, alignment="latest_align"):
        """ Convert alignment data to phylip format data type.

            i.e.
            phylip_data = [
                ('name0', 'ASDGASHASASDFASGASDFAS'),
                ('name1', 'JSERHSNDGAJEHDFDSGSHRE'),
                ...
            ]
        """
        phylip_data = []
        for id, homolog in self._Alignment._HomologSet.homologs.items():
            phylip_data.append((id, getattr(homolog, alignment)))
        return phylip_data

    def _align_to_metadata(self, tags=None):
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
        sequence_metadata = []
        if tags is None:
            for id, homolog in self._Alignment._HomologSet.homologs.items():
                sequence_metadata.append(h.attrs)
        else:
            for id, homolog in self._Alignment._HomologSet.homologs.items():
                metadata = dict([(t, getattr(homolog, t)) for t in tags])
                sequence_metadata.append(metadata)
        return sequence_metadata

    @write_to_file
    def fasta(self, alignment="latest_align", tags=("id",)):
        """ Write alignment to fasta format. """
        data = self._align_to_sequence_data(alignment=alignment, tags=tags)
        output = fasta.write(data)
        return output

    @write_to_file
    def phylip(self, alignment="latest_align"):
        """" Write alignment to phylip format. """
        data = self._align_to_phylip_data(alignment=alignment)
        output = phylip.write(data)
        return output

    @write_to_file
    def csv(self, tags=None):
        """ Write alignment to csv format. """
        sequence_metadata = self._align_to_metadata(tags=None)
        output = csv.write(sequence_metadata)
        return output

class Read(object):

    def __init__(self, Alignment):
        """ Module for reading alignment for Homologset from various formats. """
        self._Alignment = Alignment

    def _sequence_data_to_alignment(self, data, tags=["id"]):
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

        # Iterate through the sequence data
        for pair in data:
            attributes = pair[0]
            alignment = pair[1]

            id = attributes[index]

            if len(attributes) != len(tags):
                raise Exception(""" Number of attributes do not match number of tags given. """)

            homolog = getattr(self._Alignment._HomologSet, id)
            homolog.add_alignment(alignment)

        return self._Alignment

    def _phylip_data_to_alignment(self, phylip_data):
        """ Add phylip data to alignment.
        """
        # Iterate through the sequence data
        for pair in data:
            id = pair[0]
            alignment = pair[1]

            if len(attributes) != len(tags):
                raise Exception(""" Number of attributes do not match number of tags given. """)

            homolog = getattr(self._Alignment._HomologSet, id)
            homolog.add_alignment(alignment)

        return self._Alignment

    @read_from_file
    def fasta(self, data):
        """ Read alignment from fasta file. """
        sequence_data = fasta.read(data)
        self._sequence_data_to_alignment(sequence_data)
        return self._Alignment

    @read_from_file
    def phylip(self, data):
        """ Read alignment from fasta file. """
        sequence_data = phylip.read(data)
        self._phylip_data_to_alignment(sequence_data)
        return self._Alignment
