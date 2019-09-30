__description__ = \
"""
Wrap cd-hit for interfacing with phylopandas.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2019-08-16"

import phylopandas as phy

import pandas as pd
import numpy as np

import subprocess, random, string, os

def run(df,cutoff=0.9,keep_tmp=False):
    """
    Run cd-hit on a phylopandas data frame to remove similar sequences.
    Sequences are sorted from longest to shortest before running cd-hit,
    ensuring that the larger sequences always end up in the final output.

    df: phylopandas data frame to de-duplicate.
    cutoff: percent sequence identity cutoff (between 0 and 1.0) at which to
            place two sequences in the same cluster
    keep_tmp: keep temporary cd-hit files.
    """

    # Sanity check, convert to string for later
    if cutoff < 0 or cutoff > 1.0:
        err = "cutoff must be a float between 0 and 1.0 (inclusive)\n"
        raise ValueError(err)
    cutoff = "{:.3f}".format(cutoff)

    # Make a temporary data frame
    tmp_df = df.copy()

    # Sort the data frame so the longest sequence occurs first
    tmp_df["seq_length"] = [len(s) for s in tmp_df.sequence]
    tmp_df = tmp_df.sort_values(by=["seq_length"],ascending=False)

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "{}_cd-hit_in.fasta".format(tmp_file_root)
    out_root = "{}_cd-hit_out.fasta".format(tmp_file_root)

    # Write data frame to fasta file
    phy.seqio.write.to_fasta(tmp_df,id_col="uid",filename=input_file)

    # Construct cd-hit command
    cmd = ['cd-hit', "-i", input_file, "-o", out_root,"-c",cutoff]


    # Run cd-hit
    try:
        run = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
        stdoutdata, stderrdata = run.communicate()
    except FileNotFoundError:
        err = "cd-hit does not appear to be in your path\n"
        raise RuntimeError(err)

    # Make sure it returned successfully
    if run.returncode != 0:
        err = "cd-hit failed\n"
        raise RuntimeError(err)

    # Parse cd-hit output
    new_seqs = phy.read_fasta("{}".format(out_root))

    # Remove temporary files
    if not keep_tmp:
        os.remove(out_root)
        os.remove("{}.clstr".format(out_root))
        os.remove(input_file)

    # Only grab sequences with uid that made it through cdhit
    return df[df.uid.isin(list(new_seqs.id))].copy()
