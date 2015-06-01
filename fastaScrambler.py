#!/usr/bin/env python
__description__ = \
"""
Scramble the order of a set of entries in a fasta file.
"""
__usage__ = "fastaScrambler.py fasta_file"
__author__ = "Michael J. Harms"
__date__ = "110714"

import sys, re, os, random

class FastaScramblerError(Exception):
    """
    Error class.
    """

    pass

def fastaScrambler(fasta_file,collapse=True):
    """
    Scramble the order of a set of entries in a fasta file.

    If collapse == True, collapse the entire sequence onto a single line.
    """

    f = open(fasta_file,'r')
   
    out = [] 
    tmp_out = []
    line = f.readline()
    while line != "": 

        # If this is a new entry...
        if line.startswith(">"):

            # If we've been recording the previous entry, record it
            if len(tmp_out) != 0:
                if collapse == True:
                    tmp_out = [o.strip() for o in tmp_out]
                    tmp_out[0] = "%s\n" % tmp_out[0]

                out.append("%s\n" % ("".join(tmp_out)))
                tmp_out = []
                
        tmp_out.append(line)

        line = f.readline()
    
    # Print the last sequence        
    if len(tmp_out) != 0:
        if collapse == True:
            tmp_out = [o.strip() for o in tmp_out]
            tmp_out[0] = "%s\n" % tmp_out[0]

        out.append("%s\n" % ("".join(tmp_out)))

    f.close()

    random.shuffle(out) 

    return out

def main(argv=None):
    """
    Main function to parse command line.
    """   
 
    if argv == None:
        argv = sys.argv[1:]

    try:
        fasta_file = argv[0]
    except IndexError:
        err = "Incorrect number of arguments!\nUSAGE:\n\n%s\n\n" % __usage__
        raise FastaScramblerError(err)

    return fastaScrambler(fasta_file) 

if __name__ == "__main__":
    print "".join(main())
