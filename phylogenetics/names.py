# -------------------------------------------------------
# Change the names in a file to any format
# -------------------------------------------------------
from os.path import splitext

def newick_names(name):
    """ Quality controll for changing names in newick file. """
    name = name.strip()
    name = name.replace("(", "[")
    name = name.replace(")", "]")
    name = name.replace(",", "-")
    return name

def fasta_names(name):
    """ """
    return name

FORMATS = {".nwk":newick_names, ".fasta":fasta_names}

def switch(filename, homologs, current_name, new_name, format=""):
    """ Switch between names. New name can be a list
        of keys in the homolog object (might be necessary
        for unique name in newick).
    """
    # Load file and split filename for later.
    f = open(filename, "r")
    text = f.read()
    f.close()

    prefix = splitext(filename)[0]
    ext = splitext(filename)[1]

    # Construct mapping for name replace.
    mapping = homologs.get_map(current_name, new_name)

    # If new_name was a list.
    if isinstance(new_name,list):
        for m in mapping:
            mapping[m] = " - ".join(mapping[m])

    # Any quality control for the given file?
    # If so, format name accordingly
    try:
        func = FORMATS[ext]
        for m in mapping:
            name = mapping[m]
            name = func(name)
            mapping[m] = name
    except:
        pass

    # Replace the old name with this new constructed name
    for key in mapping:
        text = text.replace(key, mapping[key])

    new_file = prefix + "_names" + ext
    f = open(new_file, "w")
    f.write(text)
    f.close()
