# --------------------------------------
# Useful tools for handling fasta files.
# --------------------------------------
import pickle

def split_fasta(master_fasta):
    """
        Split a fasta into multiple single sequence fastas.
    """
    f = open(master_fasta, 'r')
    lines = f.readlines()
    f.close()
    # find the indices of all starting sequences
    indices = [i for i in range(len(lines)) if lines[i][0] == ">"]
    indices += [len(lines)]
    files = list()
    for i in range(len(indices)-1):
        filename = "tmp_seq" + str(i)
        files.append(filename)
        filename += ".fasta"
        f = open(filename, "w")
        f.write("".join(lines[indices[i]:indices[i+1]]))
    return files

def concatenate_files(filenames, output):
    """ Concatenate many files. """
    with open(output, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def load_homologset(filename):
    """ Load a homologset. """
    f = open(filename, "rb")
    homologset = pickle.load(f)
    f.close()
    return homologset

def get_fasta_names(filename):
    """ Get everthing after the `>` in a fasta file (without the sequence)"""
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    names = list()
    for l in lines:
        # Find start of sequence data
        if l[0] == ">":
            # Append this line to names list, stripping all white space
            # after the last character.
            names.append(l[1:].rstrip())

    return names
