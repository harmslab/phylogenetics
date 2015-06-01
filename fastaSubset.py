#!/usr/bin/env python
__description__ = \
"""
Take a subset of the entries in a fasta file where something in the header 
matches "pattern". -e forces an exact match to everything after >.
"""
__usage__ = "fastaSubset.py fasta_file pattern OR file_with_patterns [-e]"
__author__ = "Michael J. Harms"
__date__ = "110604"

import sys, re, os

class FastaSubsetError(Exception):
    """
    Error class.
    """

    pass

def fastaSubset(fasta_file,patterns,collapse=True,exact=False,inverse=False):
    """
    Take a subset of entries in a fasta file according to whether or not the
    header contains one of the patterns in patterns.  If collapse == True,
    collapse the entire sequence onto a single line.
    """

    compiled_patterns = [re.compile(p) for p in patterns]

    f = open(fasta_file,'r')
    
    out = []
    recording = False
    line = f.readline()
    while line != "": 

        # If this is a new entry...
        if line.startswith(">"):

            # If we've been recording the previous entry, print it out
            if len(out) != 0:
                if collapse == True:
                    out = [o.strip() for o in out]
                    out[0] = "%s\n" % out[0]
                    print "".join(out)
                else:
                    print "".join(out),

                sys.stdout.flush()

                out = []

            # Now see if this entry is matches one of the patterns
            recording = False
            if exact:
                line_search = "%s\n" % line
            else:
                line_search = "%s" % line

            if not inverse:
                for p in compiled_patterns:
                    if p.search(line_search) != None:
                        recording = True
                        break
            else:
                recording = True
                for p in compiled_patterns:
                    if p.search(line_search) != None:
                        recording = False
                        break
                
                
        # If we're recording, er, record
        if recording == True:
            out.append(line)

        line = f.readline()
    
    # Print the last sequence        
    if len(out) != 0:
        if collapse == True:
            out = [o.strip() for o in out]
            out[0] = "%s\n" % out[0]
            print "".join(out)
        else:
            print "".join(out),

    f.close()

def main(argv=None):
    """
    Main function to parse command line.
    """   
 
    if argv == None:
        argv = sys.argv[1:]

    try:
        fasta_file = argv[0]
        pattern = argv[1]
    except IndexError:
        err = "Incorrect number of arguments!\nUSAGE:\n\n%s\n\n" % __usage__
        raise FastaSubsetError(err)

    # Grab optional "exact match" argument
    try:
        if argv[2] == "-e":
            exact = True
        else:
            err = "Argument %s not recognized!\n" % argv[2]
            raise FastaSubsetError(err)
    except IndexError:
        exact = False

    if os.path.exists(pattern):
        f = open(pattern,'r')
        lines = f.readlines()
        lines = [l.split("#")[0].strip() for l in lines]
        if exact:
            patterns = ["%s\n" % l for l in lines if l != ""]
        else:
            patterns = [l for l in lines if l != ""]
        if len(patterns) == 0:
            err = "Specified pattern file (%s) does not contain any lines!\n" 
            err = "%s" % err
            raise FastaSubsetError(err)
    else:
        patterns = [pattern]
   
 
    fastaSubset(fasta_file,patterns,exact=exact) 

if __name__ == "__main__":
    main()
