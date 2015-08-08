# -------------------------------------------------------
# Change the names in a file to any format
# -------------------------------------------------------
from os.path import splitext

def switch(file, homologs, current_name, new_name):
    """ Switch between names. New name can be a list
        of keys in the homolog object (might be necessary
        for unique name in newick).    
    """
    f = open(file, "r")
    text = f.read()
    f.close()

    # Construct mapping for name replace.
    mapping = homologs.get_map(current_name, new_name)
    
    # If new_name was a list.
    if isinstance(new_name,list):
        for m in mapping:
            mapping[m] = " - ".join(mapping[m])
    
    # Replace the old name with this new constructed name
    for key in mapping:
        text = text.replace(key, mapping[key])
    
    new_file = splitext(file)[0] + "_names.fasta"
    f = open(new_file, "w")
    f.write(text)
    f.close()
    
    
