"""Module for reading and writing csv strings.
"""

def read(data, delimiter=","):
    """ Read a fasta string.

        Returns a list of tuple pairs. First value is header tags. Second
        value is sequence.
    """
    # Split data by rows
    rows = data.split("\n")

    # First row is table header
    header = rows[0]
    # Further rows are sequence data
    sequences = rows[1:]

    # Get attributes from header line
    tags = header.split(delimiter)

    # Construct metadata for each sequence
    sequence_metadata = []
    for s in sequences:
        values = s.split(delimiter)
        metadata = dict([(tags[i], values[i]) for i in range(len(tags))])
        sequence_metadata.append(metadata)

    return sequence_metadata

def write(sequence_metadata, tags=None, delimiter=","):
    """ Write a fasta string.

        Arguments:
        ---------

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
    if tags is None:
        # Get tags for each sequence
        tags = list(sequence_metadata[0].keys())

    # Header to csv string
    data = delimiter.join(tags) + "\n"

    for s in sequence_metadata:
        string = ",".join([s[t] for t in tags]) + "\n"
        data += string

    return data.strip()
