from py_blast import blast
from fasta_tools import split_fasta

def blast_seeds(input, output):
    """Blast a fasta file of seed sequences to build rough master set of sequences."""

    # Split fastas into individual files
    fastas = split_fasta(input)
    processes = [blast(f + ".fasta", f +"_blast.txt") for f in fastas]


def main():


if __name__ == "__main__":
    main()
