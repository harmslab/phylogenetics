#!/usr/bin/env python
__description__ = \
"""
Take the output from a codeml reconstruction of ancestral states and generate
ancestor(s) sampled from the posterior probabilty distribution recorded in 
that file.
"""
__author__ = "Michael J. Harms"
__date__ = "120119"
__usage__ = "sampleAncestorPosterior.py codeml_dat_file [num_ancestors = 1]"

import sys
from random import random
from math import log

def readDatFile(dat_file,take_only_top=20):
    """
    Read a codeml dat file.

    take_only_top is an integer that forces the program to only sample the top
    N states.  By default it is set to 20, so all states will be sampled.
    """

    f = open(dat_file,'r')
    lines = f.readlines()
    f.close()

    # Remove blank lines and commented (#) lines
    lines = [l for l in lines if l[0] != "#" and l.strip() != ""]

    # Read the amino acids and posterior probabilities for each site into a
    # list.
    aa_list = []
    pp_list = []
    for l in lines:
        aa_list.append([])
        pp_list.append([])
        columns = l.split()[1:]
        for i in range(0,len(columns),2):

            if take_only_top and i > (2*take_only_top):
                break
                
            aa_list[-1].append(columns[i])
            pp_list[-1].append(float(columns[i+1]))
            
        # Normalize the values in pp to be between 0 and 1
        pp_list[-1] = [p/sum(pp_list[-1]) for p in pp_list[-1]]


    # Combine amino acid and posterior probability data into a single list
    anc_data = [aa_list,pp_list]

    return anc_data

def samplePosterior(anc_data):
    """
    Sample the posterior probability data for an ancestor and generate a 
    sequence.
    """

    aa_list = anc_data[0]
    pp_list = anc_data[1]

    maxL = 1.0
    resampledL = 1.0
    diff = 0
    mean_p = 0
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

        mean_p = mean_p + pp[i] 
        maxL = pp[0]*maxL
        resampledL = pp[i]*resampledL
        if i != 0:
            diff = diff + 1
 
        sequence.append(aa[i]) 

    sequence = "".join(sequence)

    return (resampledL, maxL, mean_p/len(aa_list), diff), sequence        
    
 

def main(argv=None):
    """
    Main function.
    """

    if argv == None:
        argv = sys.argv[1:]
   
    try:
        dat_file = argv[0]
    except IndexError:
        err = "Incorrect number of arguments!\n\nUsage:\n\n%s\n\n" % __usage__
        raise IOError(err)

    try:
        num_ancestors = int(argv[1])
    except IndexError:
        num_ancestors = 1
    except ValueError:
        err = "Invalid argument (%s)!\n\nUsage:\n\n%s\n\n" % (argv[1],__usage__)
        raise IOError(err)

    anc_data = readDatFile(dat_file)

    out = []
    root = dat_file.split(".")[0]
    for i in range(num_ancestors):
        stats, sequence = samplePosterior(anc_data)
        header = ">%s %10i%10.3f%10.3f%10.3f%10i" % (root,i,
                                                     log(stats[0]),log(stats[1]),
                                                     stats[2],stats[3])
        out.append("%s\n%s\n\n" % (header,sequence))

    return out 
 
        

if __name__ == "__main__":
   
    # If invoked from command line, print output to stdout
 
    out = main()
    print "".join(out)
