#!/usr/bin/env python

# ----------------------------------------------------
# Blast a set of seed sequences in fasta file.
# ----------------------------------------------------

import argparse
from phylogenetics.names import switch

def main():
    """
        Run reverse blast.
    """
    # Build arguments and grab arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="File to change names.")
    parser.add_argument("homologs", help="homologs .pickle that holds info for sequences")
    parser.add_argument("old_name", help="Type of previous name")
    parser.add_argument("new_name", help="Type of new name to use.")
    args = parser.parse_args()

    # Run reverse blast
    switch(args[0], args[1], args[2], args[3])

if __name__ == "__main__":
    main()