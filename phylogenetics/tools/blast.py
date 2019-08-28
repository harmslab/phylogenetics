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


def ncbi_blast(df,db="nr",blast_program="blastp",
                keep_tmp=False,
                hitlist_size=100,
                entrez_query='(none)',
                e_value_cutoff=0.01,
                gapcosts=(11,1)):

    gaps = '{} {}'.format(*gapcosts)

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    out_file = "{}_blast_out.xml".format(tmp_file_root)

    uid = df.uid.tolist()
    sequences = df.sequence.tolist()
    this_query = [">{}\n{}\n".format(uid[i],sequences[i]) for i in range(len(uid))]

    # Run BLAST and download input
    result = NCBIWWW.qblast(blast_program, db, this_query,
                            hitlist_size=hitlist_size,
                            entrez_query=entrez_query,
                            expect=e_value_cutoff,
                            gapcosts=gaps)

    # Write temporary output
    f = open(out_file,"w")
    f.write(result.read())
    f.close()

    out_df = phy.read_blast_xml(out_file)

    if not keep_tmp:
        os.remove(out_file)

    return out_df



def local_blast(df,db,blast_program="blastp",
                 keep_tmp=False,
                 hitlist_size=100,
                 entrez_query='(none)',
                 e_value_cutoff=0.01,
                 gapcosts=(11,1),
                 **kwargs):
    """
    Perform a blast query using sequences in the data frame against the subject
    database.
    """

    gaps = '{} {}'.format(*gapcosts)

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
                 outfmt=5,
                 hitlist_size=hitlist_size,
                 entrez_query=entrez_query,
                 expect=e_value_cutoff,
                 gapcosts=gaps,
                 **kwargs)()

    out_df = phy.read_blast_xml(out_file)

    if not keep_tmp:
        os.remove(input_file)
        os.remove(out_file)

    return out_df

def run(df,db="nr",blast_program="blastp",
        keep_tmp=False,
        hitlist_size=100,
        entrez_query='(none)',
        e_value_cutoff=0.01,
        gapcosts=(11,1),
        **kwargs):
    """
    Perform a blast query either locally or against an NCBI server.

    ### HACK HACK HACK
    # At moment, just runs ncbi blasting.  Local blast function works great,
    # but need to figure out how to specify local vs. ncbi blast.  Maybe
    # argument?  Maybe figure it out from reserved ncbi db names?
    """

    return ncbi_blast(df,db,blast_program,
                      keep_tmp,
                      hitlist_size,entrez_query,e_value_cutoff,
                      gapcosts,**kwargs)
