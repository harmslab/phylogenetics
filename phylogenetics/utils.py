# --------------------------------------
# Useful tools for handling fasta files.
# --------------------------------------
import pickle
import subprocess

def run_subprocess(base, *args, **kwargs):
    """ Run a subprocess command with given set of args and kwargs.
        and handle errors.
    """ 
    f = [base]
    # Add positional arguments
    for a in args:
        f.append(a)
    # Add keyword arguments
    for kw in kwargs:
        f.append("-" + kw)
        f.append(kwargs[kw])

    # Run msaprobs using args.
    run = subprocess.Popen(f,stdout=subprocess.PIPE,stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
    #Answer Y to continue after initial pass. May need more quality check for this.
    run.stdin.write(bytes('Y', "ascii"))
    stdoutdata, stderrdata = run.communicate()

    # Check if alignment worked correctly.
    if run.returncode != 0:
        print(stdoutdata)
        err = base + " failed!\n"
        raise Exception(err)
    
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


def read_fasta(filename):
    """ Reads in a fasta file, returning a dict with defline 
        as key and sequence as value.
    """
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    
    keys = []
    vals = []
    for l in lines:
        val = ""
        if l[0] == ">":
            keys.append(l[1:].strip())
            vals.append("")
        else:
            vals[-1] += "".join(l.strip())
    return dict(zip(keys, vals))
