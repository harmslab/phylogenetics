# Useful functions for handling queries to NCBI Blast.
import os
import numpy as np
import glob
from subprocess import call
from collections import OrderedDict

from phylogenetics.base import Homolog, HomologSet
from phylogenetics.utils import split_fasta
from phylogenetics.formats import parse_blast_fasta, parse_blast_XML, DEFAULTS

# ----------------------------------------------------
# Python interface for standard commandline call to
# NCBI Blast servers
# ----------------------------------------------------

def ncbi(fasta_input, output, **kwargs):
    """
        Construct the NCBI Blast command for blast+ application. Send this command
        to the subprocess command.
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

def seeds(fasta, as_homologset=True, rm_blast=False, **kwargs):
    """ Blast a set of seed sequences.

        Arguments:
        ---------
        fasta : str
            filename for fasta containing seed sequences.
        as_homologset: bool [default=true]
            Convert blast results to homolog set.

        kwargs are passed to blasting method.
    """
    # grab just the name
    filename = os.path.splitext(fasta)[0]
    fastas = split_fasta(fasta)
    print("Total number of sequences before blasting: " + str(len(fastas)))

    # Make a directory for storing the blast results.
    cwd = os.getcwd()
    blastpath = os.path.join(cwd, "blast")
    os.mkdir(blastpath)

    # Run individual blasts
    outnames = []
    for f in fastas:
        # Make filenames
        iname = f+".fasta"
        oname = os.path.join(blastpath, f+"_blast.txt")

        # Send query to NCBI
        process = ncbi(iname, oname, kwargs)
        outnames.append(oname)

    # If homolog_set should be made, return homolog_set
    if as_homologset:
        # Convert to homologset
        homologset = to_homologset(outnames, tag_list=DEFAULTS)
        return homologset

def to_homologset(filenames, tag_list=DEFAULTS):
    """ Turn multiple blast hit XML files into homolog object

        Arguments:
        ---------
        filenames: string or list
            list of filenames for several blast results to compile into a homologset
        tag_list: tuple
            XML tags to strip from blast results and make attributes in homolog objects

        Returns:
        -------
        homologset: HomologSet object
    """
    # If only given one filename, don't stress! just load that file
    if type(filenames) != list:
        filenames = [filenames]

    hits = []
    for f in filenames:
        # Don't discriminate on the type of Blast XML file format, try both.
        try:
            # First, try blasting hits format.
            hits += parse_blast_XML(f, tag_list=tag_list)
        except:
            # Then, try blast fasta format
            hits += parse_blast_fasta(f)

    homologs = []
    for i in range(len(hits)):
        # Make a unique id for each sequence.
        unique_id = "XX%08d" % i

        # Create homolog instance for each sequence.
        homologs.append(Homolog(unique_id, **hits[i]))

    # Build homolog set from xml file
    homologset = HomologSet(homologs)
    return homologset

# ------------------------------------------------------
# Reverse Blasting Tools
# ------------------------------------------------------

def name_variants(name):
    """
        Create all variants of camelcase names.
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
