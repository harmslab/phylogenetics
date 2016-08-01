# Output from PAML
from __future__ import absolute_import

import re

def read(data):
    """
        Read PAML output, return a newick tree with ancestors labeled and
        ancestor_data in nested dictionary.
    """
    # Match tree with ancestor labels
    regex = re.compile("tree with node labels for Rod Page's TreeView\n.+\n")
    match = regex.search(data)
    tree = match.group().split("\n")[1]

    # Initialize the ancestor_data dictionary
    ancestor_data = {}

    # Compile a regular expression to find blocks of data for internal nodes
    node_regex = re.compile("""Prob distribution at node [0-9]+, by site[-\w():.\s]+\n\n""")
    # Strip the node number from this block of data.
    node_num_regex = re.compile("[0-9]+")

    # Iterate through each block of internal node data
    for node in node_regex.findall(data):
        # Initialize a dictionary for site data
        site_data = {}

        # Compile regex for matching site data
        site_regex = re.compile("(?:\w\(\w.\w{3}\) )+")

        site_num = 0
        # Iterate through each match for site data.
        for site in site_regex.findall(node):
            # Iterate through residues
            residue_data = dict([(site[i], float(site[i+2:i+7])) for i in range(0, len(site), 9)])
            # Add data to site_data
            site_data[site_num] = residue_data
            site_num += 1

        # Add site data to ancestor_data
        node_num = node_num_regex.search(node) # Get node_number defined by PAML
        ancestor_data[int(node_num.group())] = site_data

    return tree, ancestor_data
