import Bio.Phylo.Applications
import pandas as pd
import numpy as np
import subprocess, random, string, os
import phylopandas as ph
import pyasr
import shutil

def run(df,
        project_dir,
        id_col='uid',
        sequence_col='sequence',
        altall_cutoff=0.2,
        aaRatefile='lg',
        **kwargs):
    """
    Run ancestral sequence reconstruction powered by PAML.
    """

    reconstruction_df = pyasr.reconstruct(
        df,
        id_col=id_col,
        working_dir=project_dir,
        sequence_col=sequence_col,
        altall_cutoff=altall_cutoff,
        aaRatefile=aaRatefile,
        **kwargs
    )

    return reconstruction_df

"""def run(df,program="lazarus",keep_tmp=False,**kwargs):

    asr_functions = {"lazarus":_lazarus}

    # Figure out which alignment function to use
    try:
        asr_function = asr_functions[program]
    except KeyError:
        err = "Ancestral sequence reconstruction program '{}' not recognized.\n\n".format(program)
        err += "Should be one of:\n"
        programs = list(asr_functions.keys())
        programs.sort()

        for p in programs:
            err += "    {}\n".format(p)

        raise ValueError(err)

    input_file = "_"

    ph.seqio.write.to_phylip(df,id_col="uid",filename="_")

    asr_function(input_file,**kwargs)

    # Parse phyml output
    new_tree = ph.read_newick("__phyml_tree.txt")
    new_tree.loc[new_tree.type == 'leaf', 'uid'] = new_tree.id

    # Add to main DataFrame
    output = df.phylo.combine(new_tree, on='uid')

    if not keep_tmp:
        os.remove(input_file)

    return output"""
