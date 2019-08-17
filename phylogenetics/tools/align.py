__description__ = \
"""
Wrap software for multiple sequence alignment for interfacing with phylopandas.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2019-08-16"

import phylopandas as phy

import pandas as pd
import numpy as np

from Bio.Align.Applications import MuscleCommandline

import subprocess, random, string, os

def _align_muscle(input_file,output_file,**kwargs):
    """
    Run muscle.
    """

    # Run muscle
    try:
        cmd = MuscleCommandline(input=input_file, out=output_file,**kwargs)
        cmd_args = str(cmd).split()
        output = subprocess.run(args=cmd_args)
    except FileNotFoundError:
        err = "muscle does not appear to be in your path\n"
        raise RuntimeError(err)

    # Make sure it returned successfully
    if output.returncode != 0:
        err = "muscle failed\n"
        raise RuntimeError(err)

def _align_msaprobs(input_file,output_file,consistency=2,iterative_refinement=10):
    """
    Run msaprobs.

    input_file: input fasta file
    output_file: output fasta file
    consistency: consistency REPS use 0 <= REPS <= 5 (default: 2) passes of
                 consistency transformation
    iterative_refinement: iterative-refinement REPS use 0 <= REPS <= 1000
                          (default: 10) passes of iterative-refinement
    """

    try:
        cmd_args = ["msaprobs",
                    "--outfile",output_file,
                    "--consistency","{}".format(consistency),
                    "--iterative-refinement","{}".format(iterative_refinement),
                    input_file]
        output = subprocess.run(args=cmd_args)
    except FileNotFoundError:
        err = "msaprobs does not appear to be in your path\n"
        raise RuntimeError(err)

    # Make sure it returned successfully
    if output.returncode != 0:
        err = "msaprobs failed\n"
        raise RuntimeError(err)


def run(df,program="muscle",keep_tmp=False,**kwargs):

    alignment_functions = {"muscle":_align_muscle,
                           "msaprobs":_align_msaprobs}

    # Figure out which alignment function to use
    try:
        alignment_function = alignment_functions[program]
    except KeyError:
        err = "Alignment program '{}' not recognized.\n\n".format(program)
        err += "Should be one of:\n"
        programs = list(alignment_functions.keys())
        programs.sort()

        for p in programs:
            err += "    {}\n".format(p)

        raise ValueError(err)

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "{}_align_in.fasta".format(tmp_file_root)
    output_file = "{}_align_out.fasta".format(tmp_file_root)

    phy.seqio.write.to_fasta(df,id_col="uid",filename=input_file)

    alignment_function(input_file,output_file,**kwargs)

    # Parse muscle output
    new_seqs = phy.read_fasta("{}".format(output_file))
    new_seqs = pd.DataFrame({"uid":new_seqs.id,
                             "sequence":new_seqs.sequence})

    # Drop the sequence information from the original frame
    to_merge = df.copy()
    to_merge = to_merge.drop(labels=["sequence"],axis=1)

    # Merge sequence information back in, this time aligned.
    output = pd.merge(to_merge, new_seqs, on=['uid'], how='left')

    # Make sure no na got introduced
    if np.sum(pd.isnull(output.sequence)) > 0:
        err = "an unknown error caused sequences to be lost during alignment.\n"
        err += "temporary files {} and {} saved.\n".format(input_file,output_file)
        raise RuntimeError(err)

    # Nuke temporary files
    if not keep_tmp:
        os.remove(output_file)
        os.remove(input_file)

    return output
