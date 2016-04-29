# -------------------------------------------------------
# Change the names in a file to any format
# -------------------------------------------------------
from os.path import splitext

def newick_names(name):
    """ Quality control for changing names in newick file. """
    name = name.strip()
    name = name.replace("(", "[")
    name = name.replace(")", "]")
    name = name.replace(",", "-")
    name = name.replace("=", "-")
    name = name.replace(":", "-")
    name = name.replace(";", "-")
    return name

def fasta_names(name):
    """ """
    return name

FORMATS = {"newick":newick_names, "fasta":fasta_names}

def switch(homologset, current_name, new_names, format="newick"):
    """ Switch between names. New name can be a list
        of keys in the homolog object (might be necessary
        for unique name in newick).
    """
    tree = homologset.tree

    # Construct mapping for name replace.
    mapping = homologset.map(current_name, new_names)

    # If new_name was a list.
    if isinstance(new_names,list):
        for m in mapping:
            mapping[m] = " - ".join(mapping[m])

    # Any quality control for the given file?
    # If so, format name accordingly
    try:
        func = FORMATS[format]
        for m in mapping:
            name = mapping[m]
            name = func(name)
            mapping[m] = name
    except:
        pass

    # Replace the old name with this new constructed name
    for key in mapping:
        tree = tree.replace(key, mapping[key])

    return tree
