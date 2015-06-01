#!/usr/bin/env python
__description__ = \
"""
renumberAncestor.py

"""
__author__ = "Michael J. Harms"
__usage__ = "comapreAncestors.py ancestor_file numbering_file"
__date__ = "100726"

import sys

class renumberAncestorError(Exception):
    """
    General error class for this module.
    """
    
    pass


def renumberAncestor(ancestor_file,renumber_file):
    """
    """

    # Create a dictionary mapping old numbering to new numbering
    f = open(renumber_file,'r')
    renum = f.readlines()
    f.close()

    renum = [l.split() for l in renum
             if l.strip() != "" and not l.startswith("#")]
    renum = dict([(int(r[0]),int(r[1])) for r in renum])

    # Read through ancestor file, renumbering away 
    f = open(ancestor_file,'r')
    anc = f.readlines()
    f.close()
  
    out = [] 
    for a in anc:
        if a[0] == "#" or a.startswith("         pos") or a.strip() == "":
            out.append(a)
        else:
            current = int(a[6:12])
            out.append("%s%6i%s" % (a[:6],renum[current],a[12:]))

    return "".join(out)
    

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        ancestor_file = argv[0]
        renumber_file = argv[1]
    except IndexError:
        err = "Incorrect number of arguments!\n\n%s\n\n" % __usage__
        raise CompareAncestorError(err)

    out = renumberAncestor(ancestor_file,renumber_file)

    print out


if __name__ == "__main__":
    main() 
