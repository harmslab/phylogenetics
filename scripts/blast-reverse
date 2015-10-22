#!/usr/bin/env python

# ----------------------------------------------------
# Simple script for running a reverse blast against an 
# organism name
# ----------------------------------------------------

import argparse

def main():
    """
        Run reverse blast.
    """
    # Build arguments and grab arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input fasta filename")
    parser.add_argument("organism", help="Organism or taxa to reverse blast against.")
    args = parser.parse_args()

    # Run reverse blast
    reverse_blast(args.input, args.organism)

if __name__ == "__main__":
    main()