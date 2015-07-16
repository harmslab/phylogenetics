# --------------------------------------
# Useful tools for handling fasta files.
# --------------------------------------

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
