__author__ = "Zach Sailer"
__date__ = "7-9-15"

import os
import numpy as np
from subprocess import call
from collections import OrderedDict
import argparse

def name_variants(name):
    """
    create all variants of camelcase names.
    """
    variants = []
    variants.append(name.upper()) # all upper case
    variants.append(name.lower()) # all lower case
    variants.append(name[0].upper() + name[1:]) # capitalize first letter
    # capitalize first letter of all words
    try:
        words = name.split(" ")
        if len(words) > 1:
            camel_words = [w[0].upper() + w[1:] for w in words]
            camel_words
            variants.append(" ".join(camel_words))
    except:
        pass
    return variants

def blast_organism(fasta_input, output, taxa):
    """
        Place NCBI against an organism using blastp command line tools.
    """
    # blastp -query input_file -out output_file -db nr -entrez_query taxa -remote
    args = OrderedDict({
        "-query": fasta_input,
        "-out": output,
        "-entrez_query": taxa,
        "-db": "nr",
        "-outfmt": "6 salltitles"
    })
    # Construct blast commandline call from arguments.
    command = ["blastp"]
    for a in args:
        command += [a, args[a]]
    # Send query to Blast database using Python subprocess
    command += ["-remote"]
    return call(command)
    #return command

def split_fasta(master_fasta):
    """
        Split a fasta into multiple single sequence fastas.
    """
    f = open(master_fasta, 'r')
    lines = f.readlines()
    f.close()
    # find the indices of all starting sequences
    indices = [i for i in range(len(lines)) if lines[i][0] == ">"]
    indices += [len(lines)]
    files = list()
    for i in range(len(indices)-1):
        filename = "tmp_seq" + str(i)
        files.append(filename)
        filename += ".fasta"
        f = open(filename, "w")
        f.write("".join(lines[indices[i]:indices[i+1]]))
    return files

def organism_in_blast(name, filename):
    """
        Is organism in blast results?
    """
    # make a variants possible for name
    variants = name_variants(name)
    found = False
    for v in variants:
        f = open(filename, "r")
        # Go through each line of file
        for line in f:
            # Is name in line? If so, stop loops and return True
            if v in line:
                found = True
                break
        f.close()
        if found:
            break
    return found

def reverse_blast(master_file, organism):
    """
        Reverse blast a fasta file of multiple sequences against the database
        of an organism.
    """
    # grab just the name
    filename = os.path.splitext(master_file)[0]
    fastas = split_fasta(master_file)
    print("Total number of sequences before blasting: " + str(len(fastas)))
    # Run individual blasts
    processes = [blast_organism(f + ".fasta", f +"_blast.txt", organism) for f in fastas]

    good_fastas = []
    for f in fastas:
        found = organism_in_blast(organism, f+ "_blast.txt")
        if found:
            good_fastas.append(f)

    g = open(filename + "_reversed.fasta", 'w')
    for fasta in good_fastas:
        f = open(fasta + ".fasta", "r")
        g.write(f.read())
        f.close()
    g.close()
    print("Final number of sequences after blasting: " + str(len(good_fastas)))

    for f in fastas:
        os.remove(f + ".fasta")
        os.remove(f + "_blast.txt")
    print("Done!")


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
