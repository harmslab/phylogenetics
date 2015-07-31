#!/usr/bin/env python
__description__ = \
"""
Convert a fasta alignment file into a phylip alignment file.  Took some tweaking,
but this precise format of phylip file can be read successfully into phyml.
"""
__author__ = "Michael J. Harms"
__date__ = "100530"
__usage__ = "./fasta2phylip.py fasta"

import sys

class Fasta2PhylipError(Exception):
    """
    General error class for this module.
    """

    pass


def fasta2phylip(lines):
    """
    Take the lines from a fasta file and return in fasta format.
    """

    lines = [l for l in lines if l.strip() != ""]

    # Read fasta file
    out = []
    for l in lines:
        if l[0] == ">":
            if len(out) != 0:
                out.append("\n")
            out.append("%-10s" % l[1:11])
        else:
            out.append(l.strip())

    out = "".join(out)
    split_out = out.split("\n")

    # Quick sanity check
    lengths = dict([(len(l),i) for (i,l) in enumerate(split_out)])
    if len(lengths) > 1:
        err = "Some sequences have different lengths!\n"
        for d in lengths.keys():
            err = "%s, %i\n" % (split_out[lengths[d]][:10],d)

        raise Fasta2PhylipError(err)

    num_seq = len(split_out)
    num_columns = lengths.keys()[0]

    to_write = ["%s\n%s\n" % (o[:10],o[10:]) for o in split_out]
    to_write = "".join(to_write)

    final_out = "%i  %i\n\n%s\n" % (num_seq,num_columns-10,to_write)

    return final_out

def main(argv=None):
    """
    Parse arguments.
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        filename = argv[0]
    except IndexError:
        print __usage__
        sys.exit()

    f = open(filename)
    lines = f.readlines()
    f.close()

    return "".join(fasta2phylip(lines))

if __name__ == "__main__":
    print main()
