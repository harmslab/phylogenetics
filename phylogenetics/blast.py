# Useful functions for handling queries to NCBI Blast.
import os
import numpy as np
from subprocess import call
from collections import OrderedDict

from phylogenetics.utils import split_fasta
from phylogenetics.formats import  blast_to_homologset

# ----------------------------------------------------
# Python interface for standard commandline call to
# NCBI Blast servers
# ----------------------------------------------------

def ncbi(fasta_input, output, **kwargs):
    """
        BLAST NCBI using blastp command line tools.
    """
    # blastp -query input_file -out output_file -db nr -entrez_query taxa -remote
    args = OrderedDict({
        "-query": fasta_input,
        "-out": output,
        "-db": "nr",
        "-outfmt": "5"
    })

    # Add arguments from blast
    for kw in kwargs:
        args["-"+kw] = kwargs[kw]

    # Construct blast commandline call from arguments.
    command = ["blastp"]
    for a in args:
        command += [a, args[a]]

    # Send query to Blast database using Python subprocess
    command += ["-remote"]

    return call(command)

def seeds(fasta, to_homologset=True, **kwargs):
    """ Blast a set of seed sequences. """
    # grab just the name
    filename = os.path.splitext(fasta)[0]
    fastas = split_fasta(fasta)
    print("Total number of sequences before blasting: " + str(len(fastas)))

    # Run individual blasts
    processes = [ncbi(f + ".fasta", f +"_blast.txt", **kwargs) for f in fastas]

    # If homolog_set should be made, return homolog_set
    if to_homologset:
        # Instance of HomologSet
        full_set = HomologSet()

        for fname in glob.glob("*_blast.txt"):
            # Convert blast results to homologSet objects
            homolog_set = blast_to_homologset(filename,
                        tag_list=("Hit_id", "Hit_def", "Hit_accession",
                        "Hit_len", "Hsp_hseq"))

            # Get list of homologs and add them to the set.
            subset = homolog_set.homologs
            full_set.add_homologs(subset)

        return full_set


# ------------------------------------------------------
# Reverse Blasting Tools
# ------------------------------------------------------

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

def reverse(master_file, organism):
    """
        Reverse blast a fasta file of multiple sequences against the database
        of an organism.
    """
    # grab just the name
    filename = os.path.splitext(master_file)[0]
    fastas = split_fasta(master_file)
    print("Total number of sequences before blasting: " + str(len(fastas)))
    # Run individual blasts
    processes = [ncbi(f + ".fasta", f +"_blast.txt", entrez_query=organism) for f in fastas]

    good_fastas = []
    bad_fastas = []
    for f in fastas:
        found = organism_in_blast(organism, f+ "_blast.txt")
        if found:
            good_fastas.append(f)
        else:
            bad_fastas.append(f)

    g = open(filename + "_reversed.fasta", 'w')
    for fasta in good_fastas:
        f = open(fasta + ".fasta", "r")
        g.write(f.read())
        f.close()
    g.close()
    print("Final number of sequences after blasting: " + str(len(good_fastas)))

    for f in bad_fastas:
        os.remove(f + ".fasta")
        os.remove(f + "_blast.txt")
    print("Done!")
