__doc__ = """A module for parsing fasta strings from various database.

``FORMATS`` is a dictionary describing different fasta header formats for
different databases

Example
-------
.. code-block::

    sequences = [
        {
            "accession" : ...,
            "sequence" : ...,
            "idtype" : ...,
            "database" : ...,
        },
        {
            "accession" : ...,
            "sequence" : ...,
            "idtype" : ...,
            "database" : ...,
        },
    ]
"""

FORMATS = {
    "gb" : {
        "header" : ["dbid", "accession", "locus"],
        "database" : "GenBank"
    },
    "emb" : {
        "header" : ["dbid", "accession", "locus"],
        "database" : "EMBL Data Library"
    },
    "dbj" : {
        "header" : ["dbid", "accession", "locus"],
        "database" : "DDBJ, DNA Database of Japan"
    },
    "pir" : {
        "header" : ["dbid", "entry"],
        "database" : "NBRF PIR"
    },
    "prf" : {
        "header" : ["dbid", "name"],
        "database" : "Protein Research Foundation"
    },
    "sp" : {
        "header" : ["dbid", "accession", "entry name"],
        "database" : "SWISS-PROT"
    },
    "pdb" : {
        "header" : ["dbid", "entry", "chain"],
        "database" : "Brookhaven Protein Data Bank"
    },
    "pat" : {
        "header" : ["dbid", "country", "number"],
        "database" : "Patents"
    },
    "bbs" : {
        "header" : ["dbid", "number"],
        "database" : "Geninfo Backbone Id"
    },
    "gnl" : {
        "header" : ["dbid", "database", "identifier"],
        "database" : "General database identifier"
    },
    "ref" : {
        "header" : ["dbid", "accession", "locus"],
        "database" : "NCBI Reference Sequence"
    },
    "lcl" : {
        "header" : ["dbid", "identifier"],
        "database" : "Local Sequence Identifier"
    },
}


def parser(fasta_string):
    """Parse a fasta string and return a list of sequence dictionaries.

    The parser attempts to

    Parameters
    ----------
    fasta_string : str
        fasta string

    Returns
    -------
    sequences : list of dictionaries
        the sequence metadata pulled from fasta file.
    """
    # Get lines from fasta
    lines = string.split("\n")
    # separate headers from sequences
    mapping = {}
    for line in lines:
        # Check if current line is a header or portion of the sequence
        if line[0] == ">":
            # Set the new header
            header = lline[1:]
            # Initialize the new string
            sequence = ""
            mapping[header] = sequence
        else:
            # Append to the sequence.
            mapping[header] += line
    # separate from
    sequences = []
    n_warnings = 0
    for header, sequence in mapping:
        # attempt to identify fasta style
        try:
            # parse the header with specific parser
            header_pieces = header.split("|")
            # get dbformat from header piece
            dbformat = header_pieces[0]
            # create the metadata dict
            metadata = dict(zip(FORMATS[dbformat]["header"], header_pieces))
            # add sequence to metadata
        except:
            n_warnings += 1
            metadata = {"fasta_header" : header}
        metadata["sequence"] = sequence
        sequences.append(metadata)
    # Return number of warnings if any were given
    if n_warnings > 0:
        raise Warning("""%d warnings were raised for unidentified fasta header formats""" % n_warnings)
    return sequences
