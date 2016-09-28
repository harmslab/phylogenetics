from __future__ import absolute_import

import re

# Regular expression patter for matching in fasta.
HEADER_REGEX = re.compile("[0-9]+\s+[0-9]+")
HEADER_VAL_REGEX = re.compile("[0-9]+")

def sequential_regex(n, l):
    """ Return a regular expression for phylip file with n sequences of length l. """
    return re.compile("[A-Z0-9 ]{10}[ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-]{45}|[A-Z0-9 \t]{1,10}\n[ABCDEFGHIJKLMNOPQRSTUVWXYZ*?-]{45}")

def read(data):
    """Read a phylip string and return an Alignment metadata dictionary.
    """
    # Strip phylip header
    header = HEADER_REGEX.search(data)
    if header is None:
        raise Exception("Phylip file appears to have no header.")
    # Parse header
    header_vals = HEADER_VAL_REGEX.findall(header.group())
    n_sequences = header_vals[0] # A string with the number of sequences in file
    l_sequences = header_vals[1] # a string with length of sequences.
    # Construct the appropriate regular expression
    regex = sequential_regex(n_sequences, l_sequences)
    # Get sequences
    sequences = regex.findall(data)
    # Construct phylip data-structure
    alignment = dict(
        type="Alignment",
        module="phylogenetics.alignments"
    )
    for s in sequences:
        sequence = s[-(int(l_sequences)):]
        name = s[-len(s):-(int(l_sequences))].lstrip().rstrip()
        metadata = dict(
            id=name,
            sequence=sequence,
            type="AlignedSequence",
            module="phylogenetics.alignments"
        )
        contents.append(metadata)
    return contents

def write(phylip_data):
    """ Write a fasta string.
    Currently, this does not write phylip files in the interleaved format.
    """
    # If not already, convert phylip data to list of pairwise tuples
    if type(phylip_data) != list and len(phylip_data[0]) != 2:
        raise Exception(""" phylip_data must be a list of tuple pairs (i.e. (id, sequence))""")

    n_sequences = len(phylip_data)
    l_sequences = len(phylip_data[0][1])

    # Make header
    data = "%s    %s\n" % (n_sequences, l_sequences)

    for s in phylip_data:
        # Get the header data from tuple pair
        header = s[0][0]

        # Quality control to make sequence name is less than 10 characters.
        if len(header)> 10:
            raise Exception(""" Sequence names must be less than 10 characters. """)

        # Get sequence data
        sequence = s[1]

        # Quality control that sequence is the same length as specified by header
        if len(sequence) != l_sequences:
            raise Exception(""" Sequence named %s does not have the right length""" % header)

        # Add phylip string to data string
        data += "%s\n%s\n" %(header, sequence)

    return data.strip()
