import argparse
from phylogenetics.utils import load_homologset
from phylogenetics.cdhit import run_cdhit

def main():

    parser = argparse.ArgumentParser(description="""
    Build a tree from homologset via PhyML.
    """)

    # command line arguments
    parser.add_argument("-i", "--input", help="Input pickle file containing homologset.", type=str)
    parser.add_argument("-o", "--output", help="Output pickle filename. (No extension needed, returns .pickle)", type=str)
    parser.add_argument("-c", "--cutoff", help="Redundancy cutoff for overlapping clusters.", type=float)
    parser.add_argument("--accession", help="Accessions to make representative in cluster.", dtype=set, default=())
    parser.add_argument("--positive", help="Keywords to make representative in cluster.", dtype=set, default=())

    # Parse the arguments
    args = parser.parse_args()

    # Open homologset and run cdhit
    hs = load_homologset(args.input)
    hs2 = run_cdhit(hs, redund_cutoff=args.cutoff, args.accession, args.positive)
    hs2.write(args.output + ".pickle", format="pickle")

if __name__ == "__main__":
    main()
