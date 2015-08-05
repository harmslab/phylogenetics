# -------------------------------------------------------
# Change the names in a file to any format
# -------------------------------------------------------
from os.path import splitext

def switch(file, homologs, current_name, new_name):
    """ Switch between names."""
    f = open(file, "r")
    text = f.read()
    f.close()

    mapping = homologs.get_map(current_name, new_name)
    for key in mapping:
        text = text.replace(key, mapping[key])
    
    new_file = splitext(file)[0] + "_names.fasta"
    f = open(new_file, "w")
    f.write(text)
    f.close()