__description__ = \
"""
Wrap Gblocks for interfacing with phylopandas.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2019-08-16"


import phylopandas as phy

import pandas as pd
import numpy as np

import math, sys, subprocess, string, random, re, os

def _qual_arg(user_value,
              python_arg_name,
              gblock_arg_name,
              allowable):
    """
    Construct and sanity check a qualitative argument to
    send to gblocks.

    user_value: value to try to send to gblocks
    python_arg_name: name of python argument (for error string)
    gblock_arg_name: name of argument in gblocks
    allowable: dictionary of allowable values mapping python to
               whatever should be jammed into gblocks
    """

    if user_value in allowable.keys():
        return "-{}={}".format(gblock_arg_name,allowable[user_value])
    else:
        err = "\n\n{} '{}' not recognized\n".format(python_arg_name,
                                                user_value)
        err += "must be one of:\n"

        allowed = list(allowable)
        allowed.sort()
        for a in allowed:
            err += "    {}\n".format(a)

        raise ValueError(err)

def _quant_arg(user_value,
               python_arg_name,
               gblock_arg_name,
               min_allowed,
               max_allowed):
    """
    Construct and sanity check a quantitative argument to
    send to gblocks.

    user_value: value to try to send to gblocks
    python_arg_name: name of python argument (for error string)
    gblock_arg_name: name of argument in gblocks
    min_allowed: minimum allowable value
    max_allowed: maximum allowable value
    """

    if user_value < min_allowed or user_value > max_allowed:
        err = "\n\nvalue of {} invalid ('{}').\n".format(python_arg_name,
                                                     user_value)
        err += "must be between {} and {}.\n".format(min_allowed,
                                                     max_allowed)
        raise ValueError(err)

    return "-{}={}".format(gblock_arg_name,int(round(user_value,0)))


def run(df,
        sequence_type="protein",
        min_conserved=None,
        min_flank=None,
        max_contig_noncon=8,
        min_block_length=10,
        allowed_gap="no",
        use_sim_matrix=True,
        min_initial_block=None,
        keep_tmp=False):
    """

    Run Gblocks on a phylopandas dataframe.

    Full docs for software are here:
    http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html#Installation

    As of this writing (August 2019), Gblocks can be installed using
    bioconda (https://bioconda.github.io/).

    sequence_type: type of sequence
                   default: protein
                   allowed: ("protein","dna","codon")
                   gblock flag: -t
    min_conserved: minimum number of sequences for conserved position
                   default: num_seqs*0.5 + 1
                   allowed: (num_seqs*0.5 + 1,num_seqs)
                   gblock flag: -b1
    min_flank: minimum number of sequences for a flank position
                   default: 0.85*num_seqs
                   allowed: (min_conserved,num_seqs)
                   gblock flag: -b2
    max_contig_noncon: maximum number of contiguous non-conserved positions
                   default: 8
                   allowed: any integer > 0
                   gblock flag: -b3
    min_block_length: minimum length of a block
                   default: 10
                   allowed: any integer > 2
                   gblock flag: -b4
    allowed_gap: allowed gap positions
                   default: "no"
                   allowed: ("no","half","all")
                   gblock flag: -b5
    use_sim_matrix: use similarity matrix
                   default: True
                   allowed: True,False (True only allowed for proteins)
                   gblock flag: -b6
    min_initial_block: minimum length of an initial block
                   default: min_block_length
                   allowed: any integer > 2
                   gblock flag: -b0
    keep_tmp: keep temporary files and write out reports
                   default: False
                   allowed: False,True
    """

    tmp_file_root = "".join([random.choice(string.ascii_letters)
                             for i in range(5)])

    input_file = "{}.fasta".format(tmp_file_root)
    output_file = "{}.fasta-gb".format(tmp_file_root)
    output_summary = "{}.fasta-gb.htm".format(tmp_file_root)
    gblocks_stdout = "{}.stdout".format(tmp_file_root)
    gblocks_stderr = "{}.stderr".format(tmp_file_root)

    phy.seqio.write.to_fasta(df,id_col="uid",filename=input_file)

    cmd = ["Gblocks",input_file]

    num_seqs = len(df.uid)

    # Parse sequence type
    a = _qual_arg(sequence_type,
                  "sequence_type",
                  "t",
                  {"protein":"p","dna":"d","codon":"c"})
    cmd.append(a)

    # Parse min conserved
    min_allowed = int(math.ceil(0.5*num_seqs) + 1)
    max_allowed = num_seqs
    if min_conserved is None:
        min_conserved = min_allowed
    a = _quant_arg(min_conserved,"min_conserved","b1",
                   min_allowed,max_allowed)
    cmd.append(a)

    # Parse min flank
    min_allowed = min_conserved
    max_allowed = sys.maxsize
    if min_flank is None:
        min_flank = int(math.ceil(0.85*num_seqs))
    a = _quant_arg(min_flank,"min_flank","b2",
                   min_allowed,max_allowed)
    cmd.append(a)

    # Parse max_contig_noncon
    a = _quant_arg(max_contig_noncon,
                   "max_contig_noncon","b3",
                   1,sys.maxsize)
    cmd.append(a)

    # Parse min_block_length
    a = _quant_arg(min_block_length,
                   "min_block_length","b4",
                   2,sys.maxsize)
    cmd.append(a)

    # Parse allowed_gap
    a = _qual_arg(allowed_gap,
                  "allowed_gap",
                  "b5",
                  {"no":"n","half":"h","all":"a"})
    cmd.append(a)

    # parse use_sim_matrix
    use_sim_matrix = bool(use_sim_matrix)
    if use_sim_matrix and sequence_type != "protein":
        use_sim_matrix = False
    a = _qual_arg(use_sim_matrix,
                  "use_sim_matrix",
                  "b6",
                  {True:"y",False:"n"})
    cmd.append(a)

    # parse min_initial_block
    if min_initial_block is None:
        min_initial_block = min_block_length
    a = _quant_arg(min_initial_block,
                   "min_initial_block","b0",
                   2,sys.maxsize)
    cmd.append(a)

    if keep_tmp:
        cmd.append("-d=y")
    else:
        cmd.append("-d=n")


    # Run gblocks
    stdout = open(gblocks_stdout,"w")
    stderr = open(gblocks_stderr,"w")

    try:
        output = subprocess.run(cmd,
                                stdout=stdout,
                                stderr=stderr)
    except FileNotFoundError:
        err = "Gblocks does not appear to be in your path\n"
        raise RuntimeError(err)

    stdout.close()
    stderr.close()

    # As far as I can tell, Gblocks *always* returns 1.  At least on my
    # MacBook (Mojave 10.14.6, bioconda build.)
    if output.returncode != 1:

        err = "Gblocks failed for some reason\n"

        err += "\n\n------stdout-----\n\n"
        err += "".join(open(gblocks_stdout).read())
        err += "\n\n------stderr-----\n\n"
        err += "".join(open(gblocks_stderr).read())

        raise RuntimeError(err)

    # Parse output
    new_seqs = phy.read_fasta(output_file)
    new_seqs = pd.DataFrame({"uid":new_seqs.id,
                             "sequence":new_seqs.sequence})

    num_old_columns = len(df.loc[0].sequence)
    num_new_columns = len(new_seqs.loc[0].sequence)
    if num_new_columns == 0:
        err = "\n\nGblocks removed all columns from alignment.\n"
        err += "Trying building a better alignment, editing it manually,\n"
        err += "or altering one or more Gblocks settings to be less aggressive.\n"
        err += "See help string on this function for details.\n\n"

        raise RuntimeError(err)

    print("Gblocks reduced the number of columns from {} to {}.\n".format(num_old_columns,
                                                                          num_new_columns))

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

    # Remove temporary files
    if not keep_tmp:
        os.remove(input_file)
        os.remove(output_file)
        os.remove(output_summary)
        os.remove(gblocks_stdout)
        os.remove(gblocks_stderr)

    return output
