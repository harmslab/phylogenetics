# Module for input and output of alignment.

from .base import read_from_file, write_to_file
from .formats import fasta

class Write(object):

    def __init__(self, Alignment):
        """ Object for reading and writing alignment data from a homolog set. """
        self._Alignment = Alignment

    def _align_to_sequence_data(self, alignment="latest_align", tags=("id")):
        """ Convert alignment data to sequence data type """
        sequence_data = list()
        for h in self._Alignment.HomologSet.homologs:
            # Construct tag list
            labels = tuple()
            for t in tags:
                labels += getattr(h, t)
            # Append tuple pair to list
            sequence_data.append(labels, getattr(h, alignment))

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
        for h in self._Alignment._HomologSet.homologs:
            phylip_data.append(h.id, getattr(h, alignment))
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
            for h in self._Alignment._HomologSet.homologs:
                sequence_metadata.append(h.attrs)
        else:
            for h in self._Alignment._HomologSet.homologs:
                metadata = dict([(t, getattr(h, t)) for t in tags])
                sequence_metadata.append(metadata)
        return sequence_metadata

    @write_to_file
    def fasta(self, alignment="latest_align", tags=("id")):
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

    def csv(self, tags=None):
        """ Write alignment to csv format. """
        sequence_metadata = self._align_to_metadata(tags=None)
        output = csv.write(sequence_metadata)
        return output

class Read(object):

    def __init__(self, Alignment):
        """ Module for reading alignment for Homologset from various formats. """s
        self._Alignment = Alignment

    def _sequence_data_to_alignment(self, data):
        """ Add sequence_data to alignment """
        


    @read_from_file
    def fasta(self, data):
        """ Read alignment from fasta file. """
