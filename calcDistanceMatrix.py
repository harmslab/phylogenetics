#!/usr/bin/env python
__description__ = \
"""

"""
__author__ = "Michael J. Harms"
__usage__ = "calcDistanceMatrix.py fasta_file"
__date__ = "091221"

import sys, phyloBase

class CalcDistanceMatrixError(Exception):
    """
    General error class for this module.
    """
    
    pass


def extractAncestor(node_file,fasta_file,alternate_states=2):
    """
    """

    f = phyloBase.FastaFile(fasta_file)
    num_characters = len(f.sequences[0].sequence)

    to_take = [False for i in range(num_characters)]
    for seq in f.sequences:
        for i, character in enumerate(seq.sequence):
            if character != "-":
                to_take[i] = True

    node = readNodeFile(node_file)

    # Set up header
    out = ["%6s%6s" % (" ","pos")]
    for i in range(alternate_states):
        out.append("%6s%6s" % (("s%i" % i),("pp%i" % i)))
    out.append("%6s\n" % "total")

    # Read reconstructed characters and posterior proabilities from node list.
    # If there are fewer alternate states at a given position, append NA (the
    # missing data character from R) and a 0 posterior probability.  Sum the
    # PP we have looked at and place at the end.
    counter = 0
    for i in range(num_characters):
        if to_take[i] == True:
            output = node[i][:2*alternate_states]
            if len(output) != 2*alternate_states:
                diff = (2*alternate_states-len(output))/2
                for j in range(diff):
                    output.extend(["NA","0.000"])

            total = sum([float(output[j]) for j in range(1,2*alternate_states,2)])

            out.append("%6i%6i" % (counter,i))
            out.extend(["%6s" % o for o in output])
            out.append("%6.3f\n" % total)
            counter += 1

    return "".join(out)

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        node_file = argv[0]
        fasta_file = argv[1]
    except IndexError:
        err = "Incorrect number of arguments!\n\n%s\n\n" % __usage__
        raise ParseAncestorError(err)

    try:
        num_alternate_states = int(argv[2])
    except (IndexError,ValueError):
        num_alternate_states = 2

    ancestor = extractAncestor(node_file,fasta_file,num_alternate_states)

    print ancestor


if __name__ == "__main__":
    main() 
