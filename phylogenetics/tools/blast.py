__description__ = \
"""
Wrap BLAST for interfacing with phylopandas.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2019-08-16"

import phylopandas as phy

# Modules for blasting, etc.
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline as _local_blast

## Modules for parsing XML
#from xml.etree import ElementTree as ET

# Modules for running sundry processes
import subprocess, re, sys, os, string, random

# Set email address here for NCBI/Entrez server
Entrez.email = "PUT_YOUR_EMAIL_ADDRESS_HERE"

# PYTHON2 code to be modified
"""
def downloadSequences(accession_list,out_file,db="protein",
                      batch_download_size=50,force=False):

    Download a list of accessions in fasta/xml format.

    accession_list: list of ncbi accesion numbers
    out_file: file in which to write output in fasta/xml format
    db: database to use for accession
    batch_download_size: size of individual download packets
    force:  True/False.  Overwrite existing download file. If False, the program
            throws a notice that an old file is being used rather than re-
            downloading.


    # check for existance of out file
    if os.path.exists(out_file):
        if force:
            print "Deleting existing download file (%s)!" % out_file
            os.remove(out_file)
        else:
            print "%s already exists.  Not downloading." % out_file
            return

    print "Posting list of unique accession numbers to NCBI...",

    # Upload the list of sequences to NCBI
    to_download = ",".join([l.strip() for l in accession_list])
    post_xml = Entrez.read(Entrez.epost(db, id=to_download))
    webenv = post_xml["WebEnv"]
    query_key = post_xml["QueryKey"]

    print "DONE.\n"

    print "Downloading sequences."

    # Now download the sequences (in fasta/xml format).
    count = len(accession_list)
    out_handle = open(out_file, "w")
    for start in range(0,count,batch_download_size):
        end = min(count, start+batch_download_size)
        print "Downloading %i to %i of %i" % (start+1,end,count)
        fetch_handle = Entrez.efetch(db=db, rettype="fasta",
                                     retmode="xml",retstart=start,
                                     retmax=batch_download_size,
                                     webenv=webenv,query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)

    out_handle.close()

def seq2blastdb(fasta_set,db_name,db_type="prot",quiet=False):

    Convert a set of fasta-type sequences into a blast database.

    fasta_set: list of sequence strings in fasta format.
    db_name: database name.
    db_type: type of database to generate
    quiet: don't print status-y things


    f = open("%s.fasta" % db_name,'w')
    f.write("".join(fasta_set))
    f.close()

    if not quiet:
        print "Creating blast database %s..." % db_name,

    # Create command to run
    cmd = "makeblastdb -in %s.fasta -dbtype %s -out %s" % (db_name,db_type,db_name)
    args = shlex.split(cmd)

    # Run command
    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdout, stderr = run.communicate()

    # Make sure it exited cleanly
    if run.poll() != 0:
        print stdout
        err = "Error running makeblastdb!\n"
        raise BlastToolsError(err)

    if not quiet:
        print "DONE."

def ncbiBlast(query_fasta,out_file,db="nr",entrez_query='(none)',
              e_value_cutoff=10.0,filter=False,hitlist_size=1000,gapcosts=(11,1),
              max_query_size=10,force=False,quiet=False):

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


    # Set up proper BLAST input
    if not filter:
        filter = "none"
    else:
        filter = "on"
    gapcosts = "%i %i" % gapcosts

    if os.path.exists(out_file):
        if force:
            if not quiet:
                print "Deleting existing blast file (%s)!" % out_file
            os.remove(out_file)
        else:
            if not quiet:
                print "%s already exists.  Not performing blast." % out_file
            return

    # Go through a set of query requests no larger than max_query_size
    for counter in range(0,len(query_fasta),max_query_size):

        this_query = "".join(query_fasta[counter:counter+max_query_size])

        count = len(query_fasta[counter:counter+max_query_size])
        if not quiet:
            print "BLASTing %i sequences against the NCBI %s database" % (count,db)
            print "e_value: %.4e" % e_value_cutoff
            print "filter low complexity: %r" % filter
            print "num hits: %i" % hitlist_size
            print "gap costs: %s" % gapcosts
            print "entrez query: \'%s\'" % entrez_query

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
            print "DONE.\n"
"""


def local_blast(df,db,blast_program="blastp",keep_tmp=False,**kwargs):

    """
    Perform a blast query using sequences in the data frame against the subject
    database.
    """

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "{}_blast_in.fasta".format(tmp_file_root)
    out_file = "{}_blast_out.xml".format(tmp_file_root)

    # Write data frame to fasta file
    phy.seqio.write.to_fasta(df,id_col="uid",filename=input_file)

    _local_blast(query=input_file,
                 cmd=blast_program,
                 db=db,
                 out=out_file,
                 outfmt=5)()

    out_df = phy.read_blast_xml(out_file)

    if not keep_tmp:
        os.remove(input_file)
        os.remove(out_file)

    return out_df
