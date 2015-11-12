# Useful functions for handling queries to NCBI Blast.
import os
import glob
from subprocess import call
from collections import OrderedDict

# XML parser import
from xml.etree import ElementTree as ET

from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.utils import split_fasta, flatten_concatenated_XML

# ----------------------------------------------------
# Defaults BLAST XML parsing tags.
# ----------------------------------------------------

DEFAULTS = ("Hit_id", "Hit_def", "Hit_len","Hit_accession", "Hsp_evalue")

# ----------------------------------------------------
# Functions for talking to NCBI servers
# ----------------------------------------------------

def get_organism(defline):
    """ Get the organism in a BLAST output defline """
    return defline[defline.find("[")+1:defline.find("]")]

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


# ----------------------------------------------------
# Various Python functions for BLASTing
# ----------------------------------------------------

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

# ----------------------------------------------------
# Functions for working with Blast results
# ----------------------------------------------------

def parse_blast_xml(xml_string, tag_list=DEFAULTS):
    """
        Parse XML file of hits returned after BLASTing a sequence
        against a database.

        Args:
        ----
        filename: str
            XML filename returned from BLAST
        tag_list: tuple
            Tuple of XML tags that will be included in output data for
            each sequence.

        Returns:
        -------
        all_hits: list of dicts
            List of sequence data for all hits in XML file (with key,values
            given by tag_list).
    """

    # Fix screwed up XML if blasts were done in series...
    blast_input = flatten_concatenated_XML(xml_string,"BlastOutput_iterations")

    # Read blast file properties (in tag_list) into a list to dump out
    blast_input = ET.XML(blast_input)

    all_hits = []

    # Navigate to proper level in XML to grab hits
    for blast_output in blast_input:
        if blast_output.tag != "BlastOutput_iterations":
            continue

        for iteration in blast_output:
            if iteration.tag != "Iteration":
                continue

            for hits in iteration:
                # Catch the parent id
                if hits.tag == "Iteration_query-ID":
                    parent_id = str(hits.text).strip()

                if hits.tag != "Iteration_hits":
                    continue

                for hit in hits:

                    # Separate the tag_list out from hits and hit specs.
                    Hit_list = [t[4:] for t in tag_list if t[0:3] == "Hit"]
                    Hsp_list = [t[4:] for t in tag_list if t[0:3] == "Hsp"]

                    # Construct a dictionary of tags stripped from XML
                    properties = dict([(h.tag[4:],str(h.text).strip()) for h in hit])
                    subset = dict((k, properties[k]) for k in Hit_list)
                    all_hits.append(subset)

                    for property in hit:

                        if property.tag == "Hit_hsps":

                            for hsp in property:

                                # Construct a dictionary of tags stripped from XML
                                hsp_properties = dict([(p.tag[4:],str(p.text).strip())
                                                       for p in hsp])
                                subset = dict((k, hsp_properties[k]) for k in Hsp_list)

                                # Append inner Hsp properties to list
                                all_hits[-1].update(subset)
                                break

    return all_hits, parent_id

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
        # Open file and read.
        yup = open(f, "r")
        f = yup.read()
        yup.close()

        # Parse the HITS xml file
        data, parent = parse_blast_xml(f, tag_list=tag_list)
        hits += data

    homologs = []
    for i in range(len(hits)):
        # Make a unique id for each sequence.
        unique_id = "XX%08d" % i

        # Strip organism name from BLAST's defline argument
        try:
            organism = get_organism(hits[i]["defline"])
        except:
            organism = "NA"

        # Create homolog instance for each sequence.
        homologs.append(Homolog(unique_id, blast_query=parent, organism=organism, **hits[i]))

    # Build homolog set from xml file
    homologset = HomologSet(homologs)
    return homologset

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
