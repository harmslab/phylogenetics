#!/usr/bin/env python
__description__ = \
"""
findDiagnosticResidues.py uses the weblogo library to calculate the sequence 
entropy at sites in two alignments, then looks for diagnostic residues using a
user-specified entropy cutoff.  Each site is either unconserved (cons_none)
conserved in alignment one (cons_one), conserved in alignment two (cons_two),
conserved in the same state in both (cons_both), or diagnostic for being in 
alignment one and alignment two (diagnostic).  The aligments must have an
identical number of columns, but can have an arbitrary number or rows.  The user
may specify an entropy cutoff on the command line.

Note that gaps from the alignment are given a frequency of -1 in the final
output.
"""

__author__ = "Michael J. Harms"
__date__ = "111007"
__usage__ = "findDiagnosticResidues.py fasta_1 fasta_2 [entropy_cutoff]"

DEFAULT_ENTROPY_CUTOFF = 2.5

import sys

try:
    import weblogolib as wl
except ImportError:
    err = "\n\nPlease install the weblogo python library!\n\n"
    err = err + "\thttp://code.google.com/p/weblogo/\n\n"

    raise ImportError(err)

# Class definitions

class FindDiagnosticResiduesError(Exception):
    """
    A general error class for this module.
    """

    pass

class AlignmentColumn:
    """
    A data structure that holds the sequence entropy output for each column
    in an alignment. 
    """

    def __init__(self,column,aa_list):
        """
        Parse the output from a seqlogo calculation for a given column.
        """

        self.name = column[0].strip()
        self.entropy = float(column[21])

        self.counts = [int(i) for i in column[1:21]]
        self.aa2counts = dict([(aa_list[i],self.counts[i])
                               for i in range(len(self.counts))])

    def getMaxSet(self):
        """
        Return the amino acid(s) that are seen most often for this column 
        along with the frequency that these amino acid(s) were seen.  The list 
        of amino acids is only greater than one if there was a tie for the most
        common amino acid.
        """

        max_counts = max(self.counts)
      
        # If max_counts is zero, the columns is a gap in the alignment
        if max_counts == 0:
            return [], -1
 
        max_set = [aa for aa in self.aa2counts.keys() 
                   if self.aa2counts[aa] == max_counts]
        max_freq = float(max_counts)/sum(self.counts)

        return max_set, max_freq 


class SeqLogo:
    """
    Class that runs and parses output from weblogo for a given alignment.
    """

    def __init__(self,fasta_file):
        """
        Load a fasta file, run the calculation, and then parse the output.
        The most important output is stored in self.columns.
        """

        self.fasta_file = fasta_file

        self.runCalculation()
        self.parseOutput()

    def runCalculation(self):
        """
        Load the fasta file and run the sequence entropy calcluation.
        """   
 
        # Calculate the sequence entropy of each column in a fasta file
        f = open(self.fasta_file,'r')
        self.data = wl.LogoData.from_seqs(wl.read_seq_data(f)) 
        f.close()

    def parseOutput(self):
        """
        Parse the output from the seqlogo calculation.
        """

        # Parse the output
        data_array = str(self.data).split("\n")
        data_array = [d.split("\t") for d in data_array]
        column_data = [d for d in data_array
                       if not d[0].startswith("#") and d[0] != ""]

        # Grab the amino acid header
        size = len(column_data[0])
        header = [d for d in data_array
                  if d[0].startswith("#") and len(d) == size]
        aa = [a.strip() for a in header[0][1:21]]

        # Put each column into an instance of AlignmentColumn
        self.columns = []
        for c in column_data:
            self.columns.append(AlignmentColumn(c,aa))


def findDiagnosticResidues(fasta1,fasta2,entropy_cutoff=DEFAULT_ENTROPY_CUTOFF):
    """
    Do the sequence entropy calcluations for each alignment and then compare the
    conserved residues between them.
    """

    # Get sequence logos for each alignment
    logo1 = SeqLogo(fasta1)
    logo2 = SeqLogo(fasta2)

    # Make sure the alignments have the same number of columns
    if len(logo1.columns) != len(logo2.columns):
        err = "\n\nFasta files must have identical numbers of columns!\n\n"
        raise FindDiagnosticResiduesError(err)

    # A dictionary for counting the number of instances of a variety different
    # conservation classes.
    cons_class_counts = {"diagnostic":0,
                         "cons_both":0,
                         "cons_none":0,
                         "cons_one":0,
                         "cons_two":0}
    
    out = []   
    for i in range(len(logo1.columns)):

        # Make sure the column names match between alignments
        n1 = logo1.columns[i].name
        n2 = logo2.columns[i].name
        if n1 != n2:
            err = "Residue numbers must match between fasta files!\n\n" 
            err = "%s and %s are not the same\n\n" % (n1,n2)
            raise FindDiagnosticResiduesError(err)

        # Is the column conserved in alignment 1? 
        e1 = logo1.columns[i].entropy
        maxset1, maxset1_f = logo1.columns[i].getMaxSet()
        if e1 < entropy_cutoff:
            maxset1 = ["NA"]
        s1 = ",".join(maxset1)

        # Is the column conserved in alignment 2? 
        e2 = logo2.columns[i].entropy
        maxset2, maxset2_f = logo2.columns[i].getMaxSet()
        if e2 < entropy_cutoff:
            maxset2 = ["NA"]
        s2 = ",".join(maxset2)
  
        # Figure out whether this is unconserved, diganostic, or conserved in
        # only one lineage. 
        if maxset1 == ["NA"] and maxset2 == ["NA"]:
            cons_class = "cons_none"
        elif maxset1 != ["NA"] and maxset2 == ["NA"]:
            cons_class = "cons_one"
        elif maxset1 == ["NA"] and maxset2 != ["NA"]:
            cons_class = "cons_two"
        else:
            if maxset1 == maxset2:
                cons_class = "cons_both"
            else: 
                cons_class = "diagnostic"
        
        cons_class_counts[cons_class] += 1

        # Append to output 
        out.append("%i\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\n" % 
                   (i,n1,cons_class,e1,e2,maxset1_f,maxset2_f,s1,s2))

    # Put in header
    out.insert(0," \tpos\tcons_class\te1\te2\tf1\tf2\taa1\taa2\n")

    # Spit out summary statistics at the bottom
    out.append("# fasta1: %s\n" % fasta1)
    out.append("# fasta2: %s\n" % fasta2)
    out.append("# entropy_cutoff: %.3f\n" % entropy_cutoff)
    out.append("# diagnostic: %i\n" % cons_class_counts["diagnostic"])
    out.append("# cons_both: %i\n" % cons_class_counts["cons_both"])
    out.append("# cons_one: %i\n" % cons_class_counts["cons_one"])
    out.append("# cons_two: %i\n" % cons_class_counts["cons_two"])
    out.append("# cons_none: %i\n" % cons_class_counts["cons_none"])


    return "".join(out)
 
  
def main(argv=None):
    """
    Main function controlling this module.
    """

    if argv == None:
        argv = sys.argv[1:]

    # Parse command line
    try:
        fasta1 = argv[0]
        fasta2 = argv[1]
    except IndexError:
        err = "\n\nInsufficient number of arguments!\n\nUSAGE: %s\n" % __usage__
        raise FindDiagnosticResiduesError(err)

    # Grab optional argument
    try:
        entropy_cutoff = float(argv[2])
    except IndexError:
        entropy_cutoff = DEFAULT_ENTROPY_CUTOFF
    except:
        err = "\n\nInvalid entropy cutoff!\n\nUSAGE: %s\n" % __usage__
        raise FindDiagnosticResiduesError(err)

    print findDiagnosticResidues(fasta1,fasta2,entropy_cutoff)


# If called from the command line, run main
if __name__ == "__main__":
    main()
