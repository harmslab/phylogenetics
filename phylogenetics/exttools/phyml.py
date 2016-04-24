# Multiple sequence alignment using MSAProbs

import os, shlex, shutil, re
import subprocess

from phylogenetics.utils import read_fasta, run_subprocess

def parse_phyml_stats(data_string):
    """ Parse phyml output. """
    # Must match two patterns: `. key : value` and `- key : value`
    option_regex = re.compile("\.[\w\t ]+:.+\n|\- [\w\t ]+:.+\n")

    # Iterate through pairs
    data = {}
    for pair in option_regex.findall(data_string):
        # Split the pair
        parse = pair.split(":")
        # Each pair has either a `. ` or `- ` in front. Remove those
        key = parse[0][2:]
        # Remove whitespace
        value = parse[1].lstrip().rstrip()
        # add to data dict
        data[key] = value

    return data

def run(fname_prefix, dtype="aa", rm_tmp=True, *args, **kwargs):
    """ Simple wrapper for running PhyML within Python.

    """
    # Create a temporary fasta file from homologs as input to PhyML.
    # Build command array and run it.
    phy_fname = "%s.phy" % fname_prefix
    stats_fname = "%s.phy_phyml_stats.txt" % fname_prefix # Default name for output from phyml
    tree_fname = "%s.phy_phyml_tree.txt" % fname_prefix # default name for output tree file

    # Construct arguments for subprocess
    stuff = ("-i", phy_fname, "-d", dtype)
    stuff += args
    run_subprocess("phyml", *stuff, **kwargs)

    # Read the tree file and add to homologset.
    with open(tree_fname, "r") as f:
        tree = f.read()

    # Read the stats file.
    with open(stats_fname, "r") as f:
        stats_string = f.read()

    # Return stats.
    stats = parse_phyml_stats(stats_string)

    # Remove alignment files if you want to just keep homologs.
    if rm_tmp:
        os.remove(phy_fname)

    return tree, stats
