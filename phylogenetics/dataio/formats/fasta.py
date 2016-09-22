__doc__ = """A module for parsing fasta strings from various database.

``FORMATS`` is a dictionary describing different fasta header formats for
different databases.

The output from the parser is a list of sequence dictionaries that can be passed
directly into the HomologSet object.

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

import warnings as _warnings

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
        "database" : "UniProtKB/Swiss-Prot"
    },
    "tr" : {
        "header" : ["dbid", "accession", "entry name"],
        "database" : "UniProtKB/TrEMBL"
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


def read(fasta_string, alignment_style=False):
    """Parse a fasta string and return a list of sequence dictionaries.

    Uses defined formats from databases above.

    Parameters
    ----------
    fasta_string : str
        fasta string
    alignment_style : bool
        if True, reads header as id and sequence as alignment sequence.

    Returns
    -------
    sequences : list of dictionaries
        the sequence metadata pulled from fasta file.
    """
    # Get lines from fasta
    lines = fasta_string.strip().split("\n")
    # separate headers from sequences
    mapping = {}
    for line in lines:
        # Check if current line is a header or portion of the sequence
        if line[0] == ">":
            # Set the new header
            header = line[1:]
            # Initialize the new string
            sequence = ""
            mapping[header] = sequence
        else:
            # Append to the sequence.
            mapping[header] += line
    # separate from
    sequences = []
    # Check for alignment_style
    if alignment_style:
        return mapping
    else:
        n_warnings = 0
        for header, sequence in mapping.items():
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
            _warnings.warn("""%d warnings were raised for unidentified fasta header formats""" % n_warnings)
            return sequences

def write(metadata, alignment_style=False):
    """Return fasta string from sequence metadata.

    Parameters
    ----------
    metadata : dict
        Sequence metadata
    alignment_style : bool [default=False]
        If True, return fasta with ids as header and sequences
    """
    # Check if only one sequence or list of sequences.
    if type(metadata) != list:
        metadata = [metadata]
    fasta_string = ""
    # Simple fasta string is just id-to-sequence
    if alignment_style:
        for d in metadata:
            header = d["id"]
            line = ">" + header + "\n" + d["sequence"] + "\n"
            fasta_string += line
    else:
        for d in metadata:
            # Try getting database
            try:
                database = d["database"]
                headers = FORMATS[database]
                header = "|".join([d[h] for h in headers])
            except KeyError:
                header = d["id"]
            line = ">" + header + "\n" + d["sequence"] + "\n"
            fasta_string += line
    return fasta_string
