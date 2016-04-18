__doc__ = """

Module for reading and writing fasta strings.

"""

import re

# Regular expression patter for matching in fasta.
REGEX = re.compile(">.+\n[A-Z\-\n]+")

def read(data):
    """ Read a fasta string.

        Returns a list of tuple pairs. First value is header tags. Second
        value is sequence.
    """
    pass

def write(phylip_data, interleaved=False):
    """ Write a fasta string.

        Arguments:
        ---------

        phylip_data = [
            (10_digit_id, sequence),
            ...
        ]
    """
    # Check that phylip data is list of pairwise tuples
    if type(sequences) != list:
        sequences = list(sequences)

    data = ""

    for s in sequences:
        # Get the header data from tuple pair
        header = s[0]

        # Quality control to make sequence name is less than 10 characters.
        if len(header)> 10:
            raise Exception(""" Sequence names must be less than 10 characters. """)

        # Get sequence data
        sequence = s[1]

        # Add phylip string to data string
        data += "%s\n%s\n" %(header_str, sequence)

    return data.strip()
