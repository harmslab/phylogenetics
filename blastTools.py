#!/usr/bin/env python
__description__ = \
"""
A set of objects and functions for blasting sequences against local and
remote servers.  Also functions for filtering sequences.  Uses Biopython.  It
also requires ncbi's blast+ tools and cdhit to be installed and visible in
the path.
"""

__author__ = "Michael J. Harms"
__date__ = "110528"
__usage__ = "orthologBlast.py seed_fasta_file"
__version__ = "0.1"

# Modules for running sundry processes
import subprocess, shlex, re, sys, os, string

# Modules for blasting, etc.
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast import NCBIWWW, NCBIXML

# Modules for parsing XML
from xml.etree import ElementTree as ET

# Set email address here for NCBI/Entrez server
Entrez.email = "harms@uoregon.edu"

# Global unique name counter
homolog_counter = 0

class BlastToolsError(Exception):
    """
    General error class for this module.
    """

    pass


class Homolog:
    """
    A class that holds homologs (minimally, accession and name; maximally,
    sequence information and name of ortholog).
    """

    def __init__(self,definition,accession):
        """
        Create a new homolog.
        """

        global homolog_counter

        # Record homolog definition and accession number
        self.definition = str(definition)
        self.accession = str(accession)
        self.rank = sys.maxint

        # Assign this homolog a unique name.
        self.unique_name = "%i" % (homolog_counter)
        self.unique_name = self.unique_name.zfill(8)
        self.unique_name = "XX%s" % self.unique_name
        homolog_counter += 1

        # Initialize other variables that we don't yet know
        self.sequence = None
        self.length = None
        self.taxid = None
        self.organism = None

        self.ortholog_name = None

    def loadSequence(self,sequence,length,taxid,organism):
        """
        Load sequence data for this homolog.
        """

        self.sequence = Seq(sequence)
        self.length = length
        self.taxid = taxid
        self.organism = organism

    def formatFasta(self,seq_name=None):
        """
        Return a fasta formatted string with either "unique_name" or seq_name
        as the name.
        """

        if seq_name != None:
            name = seq_name
        else:
            seq_name = self.unique_name

        return ">%s\n%s\n" % (seq_name,self.sequence)

    def formatTabDelim(self):
        """
        Write the data for this ortholog out onto a single, tab-delimited
        line.
        """

        pretty_name = "%s-%s-%s" % (self.ortholog_name,self.organism,
                                    self.accession)

        to_write = (self.unique_name,
                    self.ortholog_name,
                    self.organism,
                    self.rank,
                    self.accession,
                    self.length,
                    str(self.sequence),
                    self.definition,
                    pretty_name)

        to_write = "\t".join(["%r" % w for w in to_write])

        return "%s\n" % to_write


def environmentCheck():
    """
    Make sure that all accessory programs are available in the path.
    """

    print("Checking for required external programs...")

    to_check = ["cdhit","blastp","tblastn","makeblastdb"]

    failed = []
    for c in to_check:
        args = shlex.split(c)
        try:
            out = subprocess.Popen(args,stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
        except OSError:
            failed.append(c)
        except:
            print("Unexpected error when running %s:" % c, sys.exc_info()[0])
            raise

    if len(failed) > 0:
        "FAILURE\n"
        return failed
    else:
        print("SUCCESS\n")

    if Entrez.email == None:
        print("No email address has been set!  To avoid typing this in the")
        print("future, edit the line 'Entrez.email = None' to point to your")
        print("email address (e.g. Entrez.email = \"A.N.Other@example.com\").")
        print("\nPlease visit the entrez website for more information about")
        print("Entrez usage rules and why providing an email address is useful.")
        print("")
        print("http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html#UserSystemRequirements")
        print("")

        email = raw_input("Please enter your email address:\n")

        Entrez.email = email

        print("")

    return []


def flattenConcatenatedXML(input_file,key_tag):
    """
    Clean up naively concatenated XML files by deleting begin/end tags that
    occur at the place where the two files were concatenated.

    NOTE: This will break and break royally if the key_tags are on the same
    lines as other important entries.
    """

    f = open(input_file,'r')
    input = f.readlines()
    f.close()

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


def parseFastaXML(sequence_file):
    """
    Load a set of sequences from an NCBI xml/fasta file into a set of Homolog
    objects.
    """

    homolog_list = []

    print("Parsing sequences in %s" % sequence_file)

    # Fix screwed up XML because sequences downloaded and output concatenated
    sequence_input = flattenConcatenatedXML(sequence_file,"TSeqSet")

    # Now we should have valid XML...
    sequence_input = ET.XML(sequence_input)
    for i, sequence in enumerate(sequence_input):
        properties = dict([(s.tag,str(s.text)) for s in sequence])

        definition = properties["TSeq_defline"].strip()
        definition = re.sub("\t"," ",definition)
        accession = properties["TSeq_gi"].strip()

        homolog_list.append(Homolog(definition,accession))
        homolog_list[-1].loadSequence(properties["TSeq_sequence"].strip(),
                                      int(properties["TSeq_length"]),
                                      int(properties["TSeq_taxid"]),
                                      properties["TSeq_orgname"].strip())

    print("DONE.")

    return homolog_list

def parseBlastXML(blast_file,tag_list=("Hit_def","Hit_id")):
    """
    Parse BLAST xml output, extracting tags specified in tag_list and putting
    into a list.  E-value is always appended after the last requested tag.
    """

    # Fix screwed up XML if blasts were done in series...
    blast_input = flattenConcatenatedXML(blast_file,"BlastOutput_iterations")

    # Read blast file properties (in tag_list) into a list to dump out
    blast_input = ET.XML(blast_input)

    all_hits = []
    for blast_output in blast_input:
        if blast_output.tag != "BlastOutput_iterations":
            continue

        for iteration in blast_output:
            if iteration.tag != "Iteration":
                continue

            for hits in iteration:
                if hits.tag != "Iteration_hits":
                    continue

                for hit in hits:
                    properties = dict([(h.tag,str(h.text)) for h in hit])
                    all_hits.append([properties[t] for t in tag_list])

                    for property in hit:
                        if property.tag == "Hit_hsps":
                            for hsp in property:
                                hsp_properties = dict([(p.tag,str(p.text))
                                                       for p in hsp])
                                all_hits[-1].append(hsp_properties["Hsp_evalue"])
                                break

    return all_hits


def downloadSequences(accession_list,out_file,db="protein",
                      batch_download_size=50,force=False):
    """
    Download a list of accessions in fasta/xml format.

    accession_list: list of ncbi accesion numbers
    out_file: file in which to write output in fasta/xml format
    db: database to use for accession
    batch_download_size: size of individual download packets
    force:  True/False.  Overwrite existing download file. If False, the program
            throws a notice that an old file is being used rather than re-
            downloading.
    """

    # check for existance of out file
    if os.path.exists(out_file):
        if force:
            print("Deleting existing download file (%s)!" % out_file)
            os.remove(out_file)
        else:
            print("%s already exists.  Not downloading." % out_file)
            return

    print("Posting list of unique accession numbers to NCBI...")

    # Upload the list of sequences to NCBI
    to_download = ",".join([l.strip() for l in accession_list])
    post_xml = Entrez.read(Entrez.epost(db, id=to_download))
    webenv = post_xml["WebEnv"]
    query_key = post_xml["QueryKey"]

    print("DONE.\n")

    print("Downloading sequences.")

    # Now download the sequences (in fasta/xml format).
    count = len(accession_list)
    out_handle = open(out_file, "w")
    for start in range(0,count,batch_download_size):
        end = min(count, start+batch_download_size)
        print("Downloading %i to %i of %i" % (start+1,end,count))
        fetch_handle = Entrez.efetch(db=db, rettype="fasta",
                                     retmode="xml",retstart=start,
                                     retmax=batch_download_size,
                                     webenv=webenv,query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)

    out_handle.close()


def runCdhit(homolog_list,redund_cutoff=0.99,tmp_file_suffix="oB_cdhit",
             keep_tmp=False):
    """
    Remove redundant homologs using cdhit.  After clustering with
    the redundancy cutoff, take the member of each cluster with the lowest
    rank.  Return a subset of homolog_list.
    """

    # Write out the fasta file with a unique name for each sequence that goes
    # >0, >1...>N.  Those numbers point to the index of the sequence in
    # homolog_list.

    # Don't do anything for empty list
    if len(homolog_list) == 0:
        print("Warning: empty list passed to cdhit!  Ignoring.")
        return homolog_list

    fasta_string = "".join([s.formatFasta(i) for i,s in enumerate(homolog_list)])
    f = open("%s.fasta" % tmp_file_suffix,'w')
    f.write(fasta_string)
    f.close()

    # Run cdhit
    cdhit_cmd = "cdhit -i %s.fasta -o %s_cdhit -c %.3f" % (tmp_file_suffix,
                                                           tmp_file_suffix,
                                                           redund_cutoff)
    args = shlex.split(cdhit_cmd)

    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdoutdata, stderrdata = run.communicate()
    if run.returncode != 0:
        print(stdoutdata)
        err = "cdhit failed!\n"
        raise BlastToolsError(err)

    # Now parse the output of cdhit and grab members of clusters with the
    # lowest rank
    f = open("%s_cdhit.clstr" % tmp_file_suffix,'r')

    out = []
    in_cluster = []
    line = f.readline()
    while line != "":

        # If we are starting a new cluster
        if line.startswith(">"):

            # ... and this is not the first cluster
            if in_cluster != []:

                # Take the member of in_cluster with the minimum rank
                ranks = [homolog_list[c].rank for c in in_cluster]
                best = in_cluster[ranks.index(min(ranks))]
                out.append(homolog_list[best])

            in_cluster = []

        # If this is not a blank line, record the seq_id in in_cluster
        elif line[0] in string.digits:
            seq_id = int(line.split(">")[1])
            in_cluster.append(seq_id)

        # Read the next line
        line = f.readline()

    # Grab the last cluster
    ranks = [homolog_list[c].rank for c in in_cluster]
    best = in_cluster[ranks.index(min(ranks))]
    out.append(homolog_list[best])

    f.close()

    # Delete temporary files
    if not keep_tmp:
        os.remove("%s.fasta" % tmp_file_suffix)
        os.remove("%s_cdhit" % tmp_file_suffix)
        os.remove("%s_cdhit.clstr" % tmp_file_suffix)
        os.remove("%s_cdhit.bak.clstr" % tmp_file_suffix)

    print("cdhit lowered redundancy @ %.3f, %i of %i kept" % (redund_cutoff,
                                                              len(out),
                                                              len(homolog_list)))

    return out

def seq2blastdb(fasta_set,db_name,db_type="prot",quiet=False):
    """
    Convert a set of fasta-type sequences into a blast database.

    fasta_set: list of sequence strings in fasta format.
    db_name: database name.
    db_type: type of database to generate
    quiet: don't print status-y things
    """

    f = open("%s.fasta" % db_name,'w')
    f.write("".join(fasta_set))
    f.close()

    if not quiet:
        print("Creating blast database %s..." % db_name)

    # Create command to run
    cmd = "makeblastdb -in %s.fasta -dbtype %s -out %s" % (db_name,db_type,db_name)
    args = shlex.split(cmd)

    # Run command
    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdout, stderr = run.communicate()

    # Make sure it exited cleanly
    if run.poll() != 0:
        print(stdout)
        err = "Error running makeblastdb!\n"
        raise BlastToolsError(err)

    if not quiet:
        print("DONE.")


def localBlast(query_fasta,out_file,db,e_value_cutoff=10.0,filter=False,
               hitlist_size=500,gapcosts=(11,1),max_query_size=10,num_threads=2,
               force=False,quiet=False):
    """
    Perform a blast using query against the subject database, generating a blast
    xml file.

    query_fasta: a list of sequences in fasta format.
    db: the BLAST database file
    out_root: the suffix for the blast output file
    e_value_cutoff: e value cutoff
    filter: apply a complexity filter mask
    hitlist_size: number of hits to grab for each input sequence
    gap_cost: tuple with (opening cost, extension cost)
    max_query_size: break the query into sets of sequences of this size and
                    combine results locally.
    force:  True/False.  Overwrite existing blast file.  If False, the program
            throws a notice that an old file is being used rather than re-
            BLASTING.
    quiet: Don't print(status-y stuff.)
    """

    tmp_in_file = "tmp_fasta_for_blast"
    tmp_out_file = "tmp_blast_output"

    if os.path.exists(out_file):
        if force:
            if not quiet:
                print("Deleting existing blast file (%s)!" % out_file)
            os.remove(out_file)
        else:
            if not quiet:
                print("%s already exists.  Not performing blast." % out_file)
            return


    # Go through a set of query requests no larger than max_query_size
    for counter in range(0,len(query_fasta),max_query_size):

        f = open(tmp_in_file,'w')
        f.write("".join(query_fasta[counter:counter+max_query_size]))
        f.close()

        count = len(query_fasta[counter:counter+max_query_size])
        if not quiet:
            print("BLASTing %i sequences against the local %s database" % (count,db))
            print("e_value: %.4e" % e_value_cutoff)
            print("filter low complexity: %r" % filter)
            print("num hits: %i" % hitlist_size)
            print("gap costs: %i %i" % gapcosts)
            print("num threads: %i" % num_threads)

        io_cmd = "blastp -query %s -db %s -outfmt 5 -out %s -num_threads %i" % \
                    (tmp_in_file,db,tmp_out_file,num_threads)
        blast_cmd = "-evalue %f -gapopen %i -gapextend %i -soft_masking %s" % \
                    (e_value_cutoff,gapcosts[0],gapcosts[1],filter)

        total_cmd = "%s %s" % (io_cmd,blast_cmd)
        args = shlex.split(total_cmd)
        run = subprocess.Popen(args,stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
        stdoutdata, stderrdata = run.communicate()
        if run.returncode != 0:
            print(stdoutdata)
            err = "blastp failed!\n"
            raise BlastToolsError(err)

        f = open(out_file, "a")
        g = open(tmp_out_file,"r")
        f.write(g.read())
        f.close()
        g.close()

        os.remove(tmp_out_file)

        if not quiet:
            print("DONE.\n")

def ncbiBlast(query_fasta,out_file,db="nr",entrez_query='(none)',
              e_value_cutoff=10.0,filter=False,hitlist_size=1000,gapcosts=(11,1),
              max_query_size=10,force=False,quiet=False):
    """
    BLAST a set of sequences against an NCBI database, generating a blast xml
    file.

    query_fasta: a list of fasta-formatted sequences to blast
    db: the BLAST database
    entrez_query: additional queries for entrez (e.g. species limit)
    e_value_cutoff: e value cutoff
    filter: apply a complexity filter mask
    hitlist_size: number of hits to grab for each input sequence
    gap_cost: tuple with (opening cost, extension cost)
    max_query_size: break the query into sets of sequences of this size and
                    combine results locally.
    force:  True/False.  Overwrite existing blast file.  If False, the program
            throws a notice that an old file is being used rather than re-
            BLASTING.
    quiet: don't print status-y things
    """

    # Set up proper BLAST input
    if not filter:
        filter = "none"
    else:
        filter = "on"
    gapcosts = "%i %i" % gapcosts

    if os.path.exists(out_file):
        if force:
            if not quiet:
                print("Deleting existing blast file (%s)!" % out_file)
            os.remove(out_file)
        else:
            if not quiet:
                print("%s already exists.  Not performing blast." % out_file)
            return

    # Go through a set of query requests no larger than max_query_size
    for counter in range(0,len(query_fasta),max_query_size):

        this_query = "".join(query_fasta[counter:counter+max_query_size])

        count = len(query_fasta[counter:counter+max_query_size])
        if not quiet:
            print("BLASTing %i sequences against the NCBI %s database" % (count,db))
            print("e_value: %.4e" % e_value_cutoff)
            print("filter low complexity: %r" % filter)
            print("num hits: %i" % hitlist_size)
            print("gap costs: %s" % gapcosts)
            print("entrez query: \'%s\'" % entrez_query)

        # Run BLAST and download input
        result = NCBIWWW.qblast("blastp", "nr", this_query,
                                hitlist_size=hitlist_size,
                                entrez_query=entrez_query,
                                expect=e_value_cutoff,
                                gapcosts=gapcosts)

        f = open(out_file, "a")
        f.write(result.read())
        f.close()
        result.close()

        if not quiet:
            print("DONE.\n")


def cleanHomologs(homolog_list,ignore=("pdb","fragment","synthetic"),
                  dubious=("putative","hypothetical","unnamed","possible",
                  "predicted","unknown","uncharacterized","mutant","isoform"),
                  gi_exclusion=(),rank_offset=0,min_length=None,max_length=None,
                  quiet=False):
    """
    Clean up sets of homologs output.  Remove duplicate, "ignore" entries,
    gi_exclusion.  Rank sequence by quality (1 by default, 2 if the hit definition
    had a word in "dubious").  You can also specify rank_offset, which allows
    you to specify a priori that this blast is better or worse than some other
    (e.g., make a nr blast 0, an EST tblastn 10).
    """

    if not quiet:
        print("Cleaning up sequences.")

    # Remove pure duplicates (with exactly the same accession number)
    homolog_list = dict([(r.accession,r) for r in homolog_list]).values()

    # Compile regular expressions
    ignore_pattern = re.compile("|".join(ignore))
    dubious_pattern = re.compile("|".join(dubious))

    gi_removal = 0
    short_removal = 0
    long_removal = 0
    ignore_removal = 0

    clean_homologs = []
    for r in homolog_list:

        # Don't keep guys in gi_exclusion
        if r.accession in gi_exclusion:
            gi_removal += 1
            continue

        # Remove short and long sequences
        if r.length != None:
            if min_length != None and r.length < min_length:
                short_removal += 1
                continue
            if max_length != None and r.length > max_length:
                long_removal += 1
                continue

        # Sanitize names (remove \t, replace with " ")
        r.definition = re.sub("\t"," ",r.definition)

        # Does one of the ignore entries occur on this line?
        tmp_definition = r.definition.lower()
        if ignore_pattern.search(tmp_definition) != None:
            ignore_removal += 1
            continue

        # Does one of the dubious entries occur on this line?
        if dubious_pattern.search(tmp_definition) != None:
            r.rank = 2 + rank_offset
        else:
            r.rank = 1 + rank_offset

        clean_homologs.append(r)

    num_rank_0 = len([h for h in clean_homologs if h.rank == 0])
    num_rank_1 = len([h for h in clean_homologs if h.rank == 1])
    num_rank_2 = len([h for h in clean_homologs if h.rank == 2])

    if not quiet:
        print("Kept %i of %i hits." % (len(clean_homologs),len(homolog_list)))
        print("gi exclusion: %i, ignored: %i, short %i, long: %i" % \
            (gi_removal,ignore_removal,short_removal,long_removal))
        print("Class 0: %i, Class 1: %i, Class 2: %i" % (num_rank_0,num_rank_1,
                                                         num_rank_2))
        print("")

    return clean_homologs


# Check the program environment when module is loaded
print("Loading blastTools version %s" % __version__)
failed_programs = environmentCheck()
if len(failed_programs) != 0:
    err = "Some required programs are not in the path!\n"
    err += "Please make sure that these program(s) are available:\n"
    err += "    \n".join(failed_programs)

    raise BlastToolsError(err)
