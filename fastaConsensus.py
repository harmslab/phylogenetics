#!/usr/bin/env python
__description__ = \
"""
fastaConsensus.py

Create a consensus sequence from a fasta file.
"""
__author__ = "Michael J. Harms"
__date__ = "2014-11-19"
__usage__ = "fastaConsensus.py fasta_file"


import sys, string

class ConsensusError(Exception):
    """
    General error class for this module.
    """

    pass

def createConsensus(fasta_file,include_gaps=True):
    """
    Create a most-common-wins consensus sequence.  If include_gaps is true, a
    '-' character can be the most common and win in the consensus.
    """

    # Read file
    f = open(fasta_file,'r')
    lines = f.readlines()
    f.close()

    # Read in all sequences
    current_seq = None
    sequences = []
    for l in lines:
        if l[0] == ">":
            if current_seq != None:
                sequences.append(list(current_seq))
            current_seq = []
        else:
            current_seq.extend(l.strip())
    sequences.append(list(current_seq))

    # Sanity check; sequences should be same length
    unique_lengths = dict([(len(s),[]) for s in sequences]).keys()
    if len(unique_lengths) != 1:
        err = "fasta file has sequences of different length!"
        raise ConsensusError(err)

    # Create list that maps indexes to letters
    index_to_aa = list(string.letters[26:])
    index_to_aa.remove("B")
    index_to_aa.remove("J")
    index_to_aa.remove("O")
    index_to_aa.remove("U")
    index_to_aa.remove("X")
    index_to_aa.remove("Z")
    if include_gaps:
        index_to_aa.append("-")

    # Create dictionary that maps letters to indexes
    aa_to_index = dict([(k,i) for i, k in enumerate(index_to_aa)])

    # Create a possible_characters x sequence_length array 
    out_array = [[0 for j in range(len(index_to_aa))]
                 for i in range(len(sequences[0]))]

    # Count how often each amino acid is seen in each column by looping over
    # seuqences
    for seq in sequences:
        for i, s in enumerate(seq):
            try:
                index = aa_to_index[s]
                out_array[i][index] += 1
            except KeyError:
                pass

    # Create the consensus sequence
    consensus = []
    for i in range(len(out_array)):
        x = [(out_array[i][j],j) for j in range(len(out_array[i]))]
        x.sort(reverse=True)

        # In this case, the first and second most common are the same.  Use the
        # first one, but report the ambiguity with a "#" at the top of the 
        # output.
        if x[0][0] == x[1][0]:
            print "# %i: %i %s, %i %s" % (i, x[0][0],index_to_aa[x[0][1]],
                                             x[1][0],index_to_aa[x[1][1]])

        consensus.append(index_to_aa[x[0][1]])

    return consensus
 
def main(argv=None):
    """
    Main function.
    """

    # Parse command line
    if argv == None:
        argv = sys.argv[1:]

    try:
        fasta_file = argv[0]
    except IndexError:
        err = "Incorrect arguments. Usage:\n\n%s\n\n" % __usage__
        raise ConsensusError(err)

    # Create consensus
    out = createConsensus(fasta_file)

    return "".join(out)

# If called from the command line, run and print to stdout
if __name__ == "__main__":
    print ">consensus"
    print(main())
