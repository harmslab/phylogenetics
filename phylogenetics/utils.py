# --------------------------------------
# Useful tools for handling fasta files.
# --------------------------------------
from __future__ import absolute_import

import os, re, pickle, subprocess, time, datetime

class SubclassError(Exception):
    """Exception raised in parent objects for methods that must be inherited."""

class LinkError(Exception):
    """Exception raised when object linking fails."""

def get_time():
    """Get a string of the date and time."""
    now = datetime.datetime.now()
    date = "%d-%d-%d_%dh%dm" % (now.month, now.day, now.year, now.hour, now.minute)
    return date

def timeit(func, *args, **kwargs):
    """ Time how long a function takes. """
    start = time.time()
    func(*args, **kwargs)
    stop = time.time()
    diff = stop-start
    return diff

def overwriting(func):
    """Wrapper function to prevent overwriting files.
    """
    def wrapper(self, fname, *args, **kwargs):
        """Check if fname exists. If it does, add number to fname.
        """
        i = 0
        path_exists = True
        while path_exists:
            path_exists = os.path.isfile(fname)
            if path_exists:
                pieces = fname.split(".")
                path = "".join(pieces[0:-1])
                ext = pieces[-1]
                new_path = path + "(" + str(i) + ")"
                fname = ".".join(new_path, ext)

        return func(self, fname=fname)

    return wrapper


def overwrite_warning(func):
    """ Wrapper function to prevent overwriting modules in HomologSet objects. """

    def wrapper(force=False, *args, **kwargs):
        """ """
        if force is True:
            return func(*args, **kwargs)
        else:
            raise Warning(
            """You are about to overwrite an object in your HomologSet object. \
            If you're okay with overwriting, call this method\ again with the \
            keyword `force=True`."""
            )

    return wrapper

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
        if len(kw) > 1:
            f.append("--"+kw)
        else:
            f.append("-" + kw)
        f.append(kwargs[kw])

    # Run command using args.
    run = subprocess.Popen(f,stdout=subprocess.PIPE,stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
    #Answer Y to continue after initial pass. May need more quality check for this.
    try:
        run.stdin.write(bytes('Y', "ascii"))
    except TypeError:
        run.stdin.write('Y')

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
