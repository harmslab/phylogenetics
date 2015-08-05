#!/usr/bin/env python

# ----------------------------------------------------
# Load xml file and create set of homologs in pickle
# ----------------------------------------------------

import argparse
import phylogenetics

def main():
    """
        Run reverse blast.
    """
    # Build arguments and grab arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input fasta filename")
    parser.add_argument("output", help="Output file (as pickle).")
    args = parser.parse_args()

    # Run reverse blast


if __name__ == "__main__":
    main()