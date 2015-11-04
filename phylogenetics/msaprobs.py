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

    return homolog_subset

def alignment_to_homologs(homolog_set, alignment_file, alignment_keys="id"):
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
    # Read the alignment file
    data = read_fasta(alignment_file)

    # Get the ids of all homologs in a set
    mapping = homolog_set.get_map(alignment_keys)
    in_homologset = list(mapping.keys())

    # Find an homologs that were removed from this homolog set in alignment
    in_alignment = list(data.keys())

    # Initialize a list of ids not in the alignment
    not_in_alignment = list()
    for h in in_homologset:
        # If this id is in homologset and not in alignment, append to list
        if h not in in_alignment:
            not_in_alignment.append(h)

    # Remove all sequences not in alignment from homolog set
    homolog_set.rm_homologs(not_in_alignment)

    ### Add alignments to each homolog in the set ###
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
        for h in homolog_set._homologs:
            h.add_attributes(** { key2 : h.latest_align })

    # Add new alignment to homolog.
    for h in homolog_set._homologs:
        h.add_attributes( **{ key:data[h.id] } )

    return homolog_set
