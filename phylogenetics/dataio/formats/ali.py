from . import fasta

def read(fasta_string):
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
    sequencelist : dict
        SequenceList metadata.
    """
    mapping = fasta.parse(fasta_string)
    # ----------------------------------------------------------------------
    # Construct a metadata object.
    # ----------------------------------------------------------------------
    sequencelist = dict(
        type="Alignment",
        module="phylogenetics.alignments",
        contents=[]
    )
    for header, sequence in mapping.items():
        # Construct the metadata.
        metadata = dict(
            header=header,
            sequence=sequence,
            module="phylogenetics.alignments",
            type="AlignedSequence"
        )
        # Add to list of sequences
        sequencelist["contents"].append(metadata)
    return sequencelist

def write(metadata):
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
    for d in metadata:
        header = d["id"]
        line = ">" + header + "\n" + d["sequence"] + "\n"
        fasta_string += line
    return fasta_string
