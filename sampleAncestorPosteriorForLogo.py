#!/usr/bin/env python
__description__ = \
"""
Resample posterior from .dat file directly, taking into account the 'interesting'
fasta file like parseAncestor.py does.  This is useful for generating sequence
logos.
"""
__author__ = "Michael J. Harms"
__usage__ = "sampleAncestorPosteriorForLogo.py node_file fasta_file [num_alternates (defaults to 1000)]"
__date__ = ""

import sys, phyloBase
from random import random

class ParseAncestorError(Exception):
    """
    General error class for this module.
    """
    
    pass

def readNodeFile(node_file):
    """
    """

    f = open(node_file,'r')
    lines = f.readlines()
    f.close()

    output = [l.split()[1:] for l in lines]

    return output 

def extractAncestor(node_file,fasta_file,alternate_states=20):
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

    # Read reconstructed characters and posterior proabilities from node list.
    # If there are fewer alternate states at a given position, append NA (the
    # missing data character from R) and a 0 posterior probability.  Sum the
    # PP we have looked at and place at the end.
    aa_list = []
    pp_list = []
    for i in range(num_characters):
        if to_take[i] == True:

            x = node[i][:2*alternate_states]
            aa_list.append([x[j] for j in range(0,len(x),2)])

            tmp_pp = [float(x[j]) for j in range(1,len(x),2)]
            pp_list.append([t/sum(tmp_pp) for t in tmp_pp])

    return aa_list, pp_list

def samplePosterior(aa_list,pp_list):
    """
    Sample the posterior probability data for an ancestor and generate a 
    sequence.
    """

    sequence = []
    for site_num in range(len(aa_list)):

        aa = aa_list[site_num]
        pp = pp_list[site_num]

        cum_sum = [pp[0]]
        for i in range(1,len(pp)):
            cum_sum.append(cum_sum[i-1]+pp[i])

        r = random()

        i = 0
        while r > cum_sum[i] and i < (len(cum_sum) - 1):
            i = i + 1

        sequence.append(aa[i]) 

    sequence = "".join(sequence)

    return sequence

def doStuff(node_file,fasta_file,num_resampled=1000,alternate_states=20):
    """
    """

    aa_list, pp_list = extractAncestor(node_file,fasta_file,alternate_states)

    out = []
    for i in range(num_resampled):
        s = samplePosterior(aa_list,pp_list)
        out.append(">seq_%i\n%s\n" % (i,s))

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
        num_resampled = int(argv[2])
    except IndexError:
        num_resampled = 1000
    
    return doStuff(node_file,fasta_file,num_resampled)


    

if __name__ == "__main__":
   
    # If invoked from command line, print output to stdout
 
    out = main()
    print "".join(out)
