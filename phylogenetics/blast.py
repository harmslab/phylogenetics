# Useful functions for handling queries to NCBI Blast.
import os
import glob
import re
from subprocess import call
from collections import OrderedDict

# XML parser import
from xml.etree import ElementTree as ET

from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.utils import split_fasta

# ----------------------------------------------------
# Defaults BLAST XML parsing tags.
# ----------------------------------------------------

DEFAULTS = ("Hit_id", "Hit_def", "Hit_len","Hit_accession", "Hsp_evalue")

# ----------------------------------------------------
# Functions for talking to NCBI servers
# ----------------------------------------------------

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
    mapping = dict([(accession_list[i], sequence_data[i]["sequence"]) for i in range(len(accession_list))])

    return mapping

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


def flatten_concatenated_XML(xml_string,key_tag):
    """
        Clean up naively concatenated XML files by deleting begin/end tags that
        occur at the place where the two files were concatenated.
        NOTE: This will break and break royally if the key_tags are on the same
        lines as other important entries.
    """
    input = xml_string.split("\n")
    set_start = re.compile("<%s>" % key_tag)
    set_end =   re.compile("</%s>" % key_tag)

    # Find all beginning and end tags...
    starts = [i for i, l in enumerate(input) if set_start.search(l) != None]

    # If this tag occurs more than once...
    if (len(starts) != 1):

        # Keep the first start reverse so we are chewing from the bottom.
        starts.pop(0)
        starts.reverse()

        # Remove all lines between each end and start, again chewing from the
        # bottom.
        for i in range(len(starts)):
            e = starts[i]
            while set_end.search(input[e]) == None:
                input.pop(e),
                e = e - 1
            input.pop(e)

    # Return freshly minted, clean XML
    return "".join(input)


def parse_blast_XML(xml_string, tag_list=DEFAULTS):
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


def parse_blast_fasta(xml_string):
    """
        Parse Blast's Fasta/XML formatted file, returning a list of each
        sequence's data.

        Args:
        ----------
        xml_string: str
            Fasta/XML formatted string from Blast Output.

        Returns:
        -------
        sequences: list of dicts
            List of sequences data in their own lists.
    """
    # Fix screwed up XML because sequences downloaded and output concatenated
    sequence_input = flatten_concatenated_XML(xml_string, "TSeqSet")

    # Now we should have valid XML...
    sequence_input = ET.XML(sequence_input)

    # XML Tag prefix to strip
    prefix = "TSeq_"

    sequences = []
    for i, sequence in enumerate(sequence_input):
        # Rip out all properties of a sequence
        properties = dict([(s.tag[len(prefix):],str(s.text).strip()) for s in sequence])

        # Append to sequences.
        sequences.append(properties)

    return sequences


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

        # Don't discriminate on the type of Blast XML file format, try both.
        try:
            # First, try blasting hits format.
            data, parent = parse_blast_XML(f, tag_list=tag_list)
            hits += data
        except Exception as e:
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
