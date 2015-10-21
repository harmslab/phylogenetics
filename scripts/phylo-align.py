import argparse
from phylogenetics.utils import load_homologset
from phylogenetics.msaprobs import run_msaprobs

def main():

    parser = argparse.ArgumentParser(description="""
    Construct a multiple sequence alignment using MSAProbs
    """)

    # command line arguments
    parser.add_argument("-i", "--input", help="Input pickle file containing homologset.", type=str)

    # Parse the arguments
    args = parser.parse_args()

    # Open homologset and run cdhit
    hs = load_homologset(args.input)
    hs2 = run_msaprobs(hs)
    hs2.write(args.output + ".pickle", format="pickle")

if __name__ == "__main__":
    main()
