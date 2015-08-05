#!/usr/bin/env python

# ----------------------------------------------------
# Blast a set of seed sequences in fasta file.
# ----------------------------------------------------

import argparse

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