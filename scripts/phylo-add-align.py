import argparse
from phylogenetics.utils import load_homologset
from phylogenetics.msaprobs import run_msaprobs

def main():

    parser = argparse.ArgumentParser(description="""
    Add an alignment to a homolog set.
    """)

    # command line arguments
    parser.add_argument("homologset", help="Input pickle file containing homologset.", type=str)
    parser.add_argument("alignment_file", help="Alignment fasta file.", type=str)

    # Parse the arguments
    args = parser.parse_args()

    # Open homologset and run cdhit
    hs = load_homologset(args.input)
    hs2 = alignment_to_homologs(hs, args.alignment_file)
    hs2.write(args.output + ".pickle", format="pickle")

if __name__ == "__main__":
    main()
