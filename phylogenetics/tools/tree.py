import Bio.Phylo.Applications
import pandas as pd
import numpy as np
import subprocess, random, string, os
import phylopandas as ph
import os

def _phyml(input_file,
           sequence_col='sequence',
           datatype='aa',
           bootstrap='-1',
           model='LG',
           frequencies='e',
           **kwargs):
    """
    Compute tree using phyml.
    """

    # Run phyml
    try:

        options = {
            'input':input_file,
            'datatype':datatype,
            'bootstrap':bootstrap,
            'model':model,
            'frequencies':frequencies,
        }

        # Update with any kwargs manually set by users.
        options.update(**kwargs)

        cmd = Bio.Phylo.Applications.PhymlCommandline(**options)
        cmd_args = str(cmd).split()
        output = subprocess.run(args=cmd_args)

    except FileNotFoundError:
        err = "phyml does not appear to be in your path\n"
        raise RuntimeError(err)

    # Make sure it returned successfully
    if output.returncode != 0:
        err = "phyml failed\n"
        raise RuntimeError(err)

def run(df, project_dir, program="phyml",keep_tmp=False,**kwargs):

    tree_functions = {"phyml":_phyml}

    # Figure out which alignment function to use
    try:
        tree_function = tree_functions[program]
    except KeyError:
        err = "Tree building program '{}' not recognized.\n\n".format(program)
        err += "Should be one of:\n"
        programs = list(tree_functions.keys())
        programs.sort()

        for p in programs:
            err += "    {}\n".format(p)

        raise ValueError(err)

    input_file = "_"

    os.chdir(project_dir)

    ph.seqio.write.to_phylip(df,id_col="uid",filename="_")

    tree_function(input_file,**kwargs)

    # Parse phyml output
    new_tree = ph.read_newick("__phyml_tree.txt")
    new_tree.loc[new_tree.type == 'leaf', 'uid'] = new_tree.id

    # Add to main DataFrame
    output = df.phylo.combine(new_tree, on='uid')

    if not keep_tmp:
        os.remove(input_file)

        os.chdir("..")

    return output
