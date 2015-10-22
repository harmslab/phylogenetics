# Useful functions for handling queries to NCBI Blast.
import os
import numpy as np
import glob
from subprocess import call
from collections import OrderedDict

from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.utils import split_fasta
from phylogenetics.formats import parse_blast_fasta, parse_blast_XML, DEFAULTS

def download(accession_list, email, out_file, db="protein",
                      batch_download_size=50, write=False):
    """
    Download the Entrez sequences from a list of accessions.

    Arguments :
    ---------
    accession_list: list
        list of ncbi accesion numbers
    out_file:
        file in which to write output in fasta/xml format
    db:
        database to use for accession
    batch_download_size:
        size of individual download packets
    force:  True/False.
        Overwrite existing download file. If False, the program
        throws a notice that an old file is being used rather than re-
        downloading.
    """
    # Biopython is imported here...  I realize this is a bit overkill for now.
    from Bio import Entrez

    Entrez.email = email

    #first get GI for query accesions
    query  = " ".join(accession_list)
    handle = Entrez.esearch( db=db,term=query,retmax=10**9 )
    giList = Entrez.read(handle)['IdList']

    #post GID list to NCBI
    search_handle = Entrez.epost(db=db, id=",".join(giList))
    search_results = Entrez.read(search_handle)
    webenv,query_key = search_results["WebEnv"], search_results["QueryKey"]

    # Now download the sequences (in fasta/xml format).
    count = len(accession_list)
    total_xml = ""
    for start in range(0,count,batch_download_size):
        end = min(count, start+batch_download_size)
        fetch_handle = Entrez.efetch(db=db, rettype="fasta",
                                     retmode="xml",retstart=start,
                                     retmax=batch_download_size,
                                     webenv=webenv,query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        total_xml += data + "\n"

    # Write to file is interested.
    if write:
        out_handle = open(out_file, "w")
        out_handle.write(data)
        out_handle.close()

    # parse the xml and get sequence data as list
    sequence_data = parse_blast_fasta(total_xml)

    # Map accession list to their sequences
    mapping = dict(zip(accession_list, sequence_data))

    return mapping


# ----------------------------------------------------
# Python interface for standard commandline call to
# NCBI Blast servers
# ----------------------------------------------------

def query(fasta_input, output, remote=True, **kwargs):
    """ Construct Blast+ application command. Send this command to subprocess.
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
    if remote:
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
        process = query(iname, oname, kwargs)
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
            data, parent = parse_blast_XML(f, tag_list=tag_list)
            hits += data
        except:
            # Then, try blast fasta format
            hits += parse_blast_fasta(f)
            parent=None

    homologs = []
    for i in range(len(hits)):
        # Make a unique id for each sequence.
        unique_id = "XX%08d" % i

        # Create homolog instance for each sequence.
        homologs.append(Homolog(unique_id, blast_query=parent, **hits[i]))

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
    processes = [query(f + ".fasta", f +"_blast.txt", entrez_query=organism) for f in fastas]

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
