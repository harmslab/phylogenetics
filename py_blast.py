# ---------------------------------------
# BLAST against NCBI database from python
# ---------------------------------------

def blast(fasta_input, output, db="nr", **kwargs):
    """
        Place NCBI against an organism using blastp command line tools.
    """
    # Using subprocess to make command at commandline
    # args constructs input command --
    # key is commandline flag, value is commandline argument
    args = OrderedDict({
        "-query": fasta_input,
        "-out": output,
        "-db": "nr",
    })
    # Add extra commands
    for k in kwargs:
        args[k] = kwargs[k]

    # Construct blast commandline call from arguments.
    command = ["blastp"]
    for a in args:
        command += [a, args[a]]
    # Send query to Blast database using Python subprocess
    command += ["-remote"]
    return call(command)

def blastp2fasta():
    """"""
    pass
