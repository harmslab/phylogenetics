__description__ = \
"""
phyloBase.py

A set of class definitions for parsing and manipulating a variety of formats
relevant to phylogenetic analysis.
"""
__author__ = "Michael J. Harms"
__date__ = "091219"
__usage__ = "not invoked from the command line"

import string

class PhyloBaseError(Exception): 
    """
    General error class for this module.
    """

    pass

class Sequence:
    """
    A class to hold sequences.
    """

    def __init__(self,sequence,header="",strict=True):
        """
        Create instance of sequence class.
        """
   
        # Record sequence and header, making sure that leading/trailing 
        # whitespace is removed. 
        self.header = header.strip()
        self.sequence = sequence.strip().upper()

        if strict:
            self.sanityCheck()

    def sanityCheck(self):
        """
        Make sure that the sequence is allowed.
        """

        legal_characters = list(string.ascii_uppercase)
        legal_characters.extend(["-","?","*"])

        not_sane = [s for s in self.sequence if s not in legal_characters]

        if len(not_sane) > 0:
            prob = ", ".join(dict([(k,[]) for k in not_sane]))
            err = "Sequence %s contains illegal characters!\n\n%s\n\n" % \
                (header,prob)
            raise PhyloBaseError(err)


class Alignment:
    """
    Class to hold an alignment (as a set of Sequence objects).
    """

    def __init__(self,sequences):
        """
        """

        self.sequences = sequences[:]
    
        # Make sure that everyone has the same length by using length as 
        # dictionary keys.  Since dictionary keys must be unique, this will lead
        # to a length-1 dictionary if all sequences have the same length.
        lengths = [len(s) for s in sequences]
        if len(dict([(l,[]) for l in lenghts])) != 1:
            err = "All sequences do not have the same length!  Bad alignment!\n"
            raise PhyloBaseError(err)

        self.length = lengths[0]


class NewickTree:
    """
    Read and write Newick-style phylogenetic trees.
    """

    def __init__(self,tree_file=None,tree_number=0):
        """
        Create an instance of NewickTree.
        """
  
        if tree_file != None:
            self.loadFile(tree_file,tree_number)
 

    def loadFile(tree_file,tree_number=0):
        """
        Read a tree file.
        """ 
        f = open(tree_file,'r')
        lines = f.readlines()
        f.close()


class FastaFile:
    """
    Read and write fasta files.
    """

    def __init__(self,fasta_file=None):
        """
        Create an instance of FastaFile.
        """           

        if fasta_file != None:
            self.loadFile(fasta_file)

    def loadFile(self,fasta_file):
        """
        Load a fasta file into memory.
        """

        f = open(fasta_file,'r')
        lines = f.readlines()
        f.close()

        # Remove trailing spaces/line endings and skip blank lines
        lines = [l.strip() for l in lines if l.strip() != ""]

        # Find all ">" entries
        entries = [i for i, l in enumerate(lines) if l.startswith(">")]
        num_seq = len(entries)

        # Grab the header names, tossing ">" and stripping spaces
        headers = [lines[entries[i]][1:].strip() for i in range(num_seq)]
       
        # Grab the actual sequences for each entry
        entries.append(-1)
        sequences = ["".join(lines[entries[i]+1:entries[i+1]])
                     for i in range(num_seq)] 
        self.sequences = [Sequence(sequences[i],headers[i])
                          for i in range(num_seq)]


