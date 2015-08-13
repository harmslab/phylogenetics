# Multiple sequence alignment using MSAProbs

import os, shlex
import subprocess

from phylogenetics.utils import read_fasta

def run_msaprobs(homolog_set, tmp_file_suffix="alignment", rm_tmp=True):
    """ Use MSAProbs to run a multiple sequence alignment on a set of homolog objects."""

    # Create a temporary fasta file from homologs as input to CDHIT.
    fname = "pre_%s.fasta" % tmp_file_suffix
    homolog_set.write(fname, format="fasta", tags=["id"])
    
    # Build msaprobs command
    msaprobs_cmd = "msaprobs pre_%s.fasta -o %s.fasta" % (tmp_file_suffix,
                                                           tmp_file_suffix)
    # Format the command into list.
    args = shlex.split(msaprobs_cmd)

    # Run msaprobs using args.
    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdoutdata, stderrdata = run.communicate()

    # Check if alignment worked correctly.
    if run.returncode != 0:
        print(stdoutdata)
        err = "msaprobs failed!\n"
        raise BlastToolsError(err)
    
    # Read the alignment file.
    data = read_fasta("%s.fasta" %tmp_file_suffix) 
    
    # Get a map of homolog ids to their object
    homolog_map = homolog_set.get_map("id")
    
    # Add sequence alignment to homolog objects in set.
    for d in data:
        homolog_map[d].add_attributes(alignment=data[d])
    
    # Remove alignment files if you want to just keep homologs.
    if rm_tmp:
        os.remove("pre_%s.fasta" % tmp_file_suffix)
        os.remove("%s.fasta" % tmp_file_suffix)

    return homolog_set

def alignment_to_homologs(homolog_set, alignment_file):
    """ Load an alignment fasta file into a set of homologs. 
    
        This function is predominantly useful for loading an alignment
        that was manually edited.
    """
    
    data = read_fasta(aligment_file)
  
    # Get a map of homolog ids to their object
    homolog_map = homolog_set.get_map("id")
    
    key = "alignment"
    counter += 0
    while hasattr(homolog_set.homologs[0], key) is True:
        key = "alignment%d" % counter
    
    for d in data:
        homolog_map[d].add_attributes({key:data[d]})

    return homolog_set