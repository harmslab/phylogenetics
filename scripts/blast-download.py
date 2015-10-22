#!/usr/bin/env python

# ----------------------------------------------------
# Blast a set of seed sequences in fasta file.
# ----------------------------------------------------

import os
import argparse
import pickle
from phylogenetics.blast import download
from phylogenetics.utils import load_homologset

def main():
    """ Download the full sequences for a set of homologs from NCBI.
    """

    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        descriptions="Downloads the full sequences for a set of homologs from NCBI.",
    )

    # Type of arguments
    parser.add_argument("input",
                        help="pickle filename containing a homolog set to download."
                        type=str)

    parser.add_argument("output",
                        help="filename to save the resulting homologset with sequences (no extension).",
                        type=str)

    parser.add_argument("email"
                        help="Email to send to Entrez. Required by NCBI's Entrez."
                        default=str)

    parser.add_argument("--db",
                        help="database to pull sequences from. See NCBI for this. [default=protein]"
                        type=str,
                        default="protein")

    parser.add_argument("--xml"
                        help="If True, will write the downloaded xml out to a file.",
                        type=bool,
                        default=False)

    # Parse the argument
    args = parser.parse_args()

    # Get the accession numbers from homologset
    hs = load_homologset(args.input)

    # Get accession ids
    accessions = [h.accession for h in hs.homologs]

    # download sequences
    sequence_data = download(accessions, args.email, args.output, db=args.db,
                          batch_download_size=50, write=args.xml)

    # Add sequences to homologset
    for h in hs.homologs:
        h.add_attributes(seq=sequence_data[h.accession])

    hs.write(args.output + ".pickle", format="pickle")

if __name__ == "__main__":
    main()
