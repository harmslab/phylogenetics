# Multiple sequence alignment using MSAProbs

import os, shlex
import subprocess

# Local imports
from phylogenetics.utils import read_fasta

def run_msaprobs(homolog_set, tmp_file_suffix="alignment", rm_tmp=True):
    """ Use MSAProbs to run a multiple sequence alignment on a set of homolog objects.
    
        Args:
        ----
        homolog_set: HomologSet object
            Instance of HomologSet class.
        tmp_file_suffix: str
            filename for temporary files
        rm_tmp: bool
            If true, removes the temporary files that it creates.
        
        Returns:
        -------
        homolog_set: HomologSet object
    """

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
        raise Exception(err)
    
    homolog_subset = alignment_to_homologs(homolog_set, "%s.fasta" %tmp_file_suffix)
    
    # Remove alignment files if you want to just keep homologs.
    if rm_tmp:
        os.remove("pre_%s.fasta" % tmp_file_suffix)
        os.remove("%s.fasta" % tmp_file_suffix)

    return homolog_set

def alignment_to_homologs(homolog_set, alignment_file):
    """ Load an alignment fasta file and add as `latest_align` attribute 
        to a set of homologs. 
    
        When doing multiple rounds of alignment, this method looks for 
        `lastest_align` attribute in homolog, replaces it, and move the old 
        alignment to a new attribute `align`+#. 
        
        Args:
        ----
        homolog_set: HomologSet object
  
        alignment_file: str
            name of aligned fasta file.
        
        Returns:
        -------
        homolog_set: HomologSet object
        
    """
    
    data = read_fasta(alignment_file)
  
    # Get a map of homolog ids to their object
    homolog_map = homolog_set.get_map("id")
    
    key = "latest_align"
    
    # If an alignment already exists in homologs, store it under
    # new attribute
    if hasattr(homolog_set.homologs[0], key): 
        # Find the last alignment
        counter = 0
        key2 = "align0"
        while hasattr(homolog_set.homologs[0], key2) is True:
            key2 = "align%d" % counter
            counter += 1
        
        # Move old alignment to new attribute in homolog
        for h in homolog_set.homologs:
            h.add_attributes(** { key2 : h.latest_align })
     
    # Add new alignment to homolog.
    for d in data:
        homolog_map[d].add_attributes( **{key:data[d]})

    return homolog_set