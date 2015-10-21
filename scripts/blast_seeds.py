#!/usr/bin/env python

# ----------------------------------------------------
# Blast a set of seed sequences in fasta file.
# ----------------------------------------------------

import os
import argparse
import pickle
from phylogenetics.blast import seeds

def main():
    """
        Blast a set of set of seed sequences from a fasta file.
    """

    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        descriptions="Blast a set of inition sequences from an input file.",
    )

    # Type of arguments
    parser.add_argument("-i", "--input",
                        help="Input fasta filename")

    parser.add_argument("-p, --pickle",
                        help="Save the blast results to a pickled HomologSet object")

    parser.add_argument("--rmblast",
                        help="Remove the blast results directory. Only works if pickle object exists."
                        default=False)

    # Parse the argument
    args = parser.parse_args()

    # If pickle option is given, set homologset to true.
    if args.pickle is not None:
        as_hs = True
    else:
        as_hs = False

    # Run reverse blast
    hs = seeds(args.input, as_homologset=as_hs, rm_blast=args.rmblast)

    # If pickle name is given, save to pickle file. 
    if as_hs:
        cwd = os.getcwd()
        fname = os.path.join(cwd, arg.pickle)
        f = open(fname, "wb")
        pickle.dump(f)
        f.close()

if __name__ == "__main__":
    main()
