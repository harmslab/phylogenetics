__author__ = "Zach Sailer"
__date__ = "7-9-15"

import os
import numpy as np
from subprocess import call
from collections import OrderedDict
import argparse
from fasta_tools import split_fasta

def blast_ncbi(fasta_input, output):
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
    # Construct blast commandline call from arguments.
    command = ["blastp"]
    for a in args:
        command += [a, args[a]]
    # Send query to Blast database using Python subprocess
    command += ["-remote"]
    return call(command)
    #return command

def split_fasta(master_fasta):
    """
        Split a fasta into multiple single sequence fastas.
    """
    f = open(master_fasta, 'r')
    lines = f.readlines()
    f.close()
    # find the indices of all starting sequences
    indices = [i for i in range(len(lines)) if lines[i][0] == ">"]
    indices += [len(lines)]
    files = list()
    for i in range(len(indices)-1):
        filename = "tmp_seq" + str(i)
        files.append(filename)
        filename += ".fasta"
        f = open(filename, "w")
        f.write("".join(lines[indices[i]:indices[i+1]]))
    return files

def full_blast(master_file):
    """
        Reverse blast a fasta file of multiple sequences against the database nr.
    """
    # grab just the name
    filename = os.path.splitext(master_file)[0]
    fastas = split_fasta(master_file)
    print("Total number of sequences before blasting: " + str(len(fastas)))
    # Run individual blasts
    processes = [blast_ncbi(f + ".fasta", f +"_blast.txt") for f in fastas]

def main():
    """
        Run reverse blast.
    """
    # Build arguments and grab arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input fasta filename")
    args = parser.parse_args()

    # Run reverse blast
    full_blast(args.input)


if __name__ == "__main__":
    main()
