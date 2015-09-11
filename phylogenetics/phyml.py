# Multiple sequence alignment using MSAProbs

import os, shlex, shutil
import subprocess

from phylogenetics.utils import read_fasta, run_subprocess

def run_phyml(homolog_set, outfile_prefix, dtype="aa", rm_tmp=False, *args, **kwargs):
    """ Construct a maximum likelihood tree using PhyML.

		__Arguments__:

		`homolog_set` : Homologset Object with homologs to construct tree.

    """
    # Create a temporary fasta file from homologs as input to CDHIT.
    fname = "%s.phy" % outfile_prefix
    homolog_set.write(fname, format="phylip")
    
    # Build command array and run it.
    stuff = ("-i", fname, "-d", dtype)
    stuff += args
    run_subprocess("phyml", *stuff, **kwargs)
    
    # Read the tree file and add to homologset.
    f = open("%s.nwk" % outfile_prefix, "r")
    tree = f.read()
    f.close()
    
    homolog_set.tree = tree
    
    # Change extension of tree to nwk
    shutil.move("%s.phy_phyml_tree.txt" % outfile_prefix, 
                "%s.phy_phyml_tree.nwk" % outfile_prefix)

    # Remove alignment files if you want to just keep homologs.
    if rm_tmp:
        os.remove("%s.phy" % outfile_prefix)

    return homolog_set
