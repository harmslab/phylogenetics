# Useful functions for handling queries to NCBI Blast.

import numpy as np

def ncbi(fasta_input, output, **kwargs):
    """
        Place NCBI against an organism using blastp command line tools.
    """
    # blastp -query input_file -out output_file -db nr -entrez_query taxa -remote
    args = OrderedDict({
        "-query": fasta_input,
        "-out": output,
        "-db": "nr",
        "-outfmt": "5"
    })

    # Add arguments from blast
    for kw in kwargs:
        args[kw] = kwargs[kw]

    # Construct blast commandline call from arguments.
    command = ["blastp"]
    for a in args:
        command += [a, args[a]]

    # Send query to Blast database using Python subprocess
    command += ["-remote"]

    return call(command)
