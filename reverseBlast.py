#!/usr/bin/env python
__description__ = \
"""
Do a primary blast with a set of seed sequences.  Take the result and then
reverse blast each entry against the entire result.  Record the top set of 
hits.  
"""

__author__ = "Michael J. Harms"
__date__ = "110528"
__usage__ = "reverseBlast.py seed_fasta_file"

import sys
from blastTools import *


def reverseBlast(seed_fasta,min_length=50,max_length=300,num_homolog_cutoff=3):
    """
    """
    
    # Read query fasta file into a set of sequences
    seed_seq = [s for s in SeqIO.parse(open(seed_fasta), format="fasta")]
    seed_list = []
    for s in seed_seq:
        seed_list.append(Homolog(s.description,None))
        seed_list[-1].rank = 0
        seed_list[-1].loadSequence(str(s.seq),len(s.seq),None,"seed")
        seed_list[-1].ortholog_name = s.description
    
    # Find maximum sequence redundancy cutoff that allows us to differentiate 
    # the seed sequence.  
    print "Optimizing false positive redundancy cutoff."
    cutoff = 0.60
    redund = []
    while len(redund) != len(seed_list):
        cutoff += 0.05
        print "Cutoff: %.2f" % cutoff
        try:
            redund = runCdhit(seed_list,redund_cutoff=cutoff)
        except BlastToolsError:
            pass 

        if cutoff > 0.95:
            err = "\nSome seed sequences have redundancy > 0.95!\n"
            err += "Run cdhit to find sequences that are too similar.\n"
            raise BlastToolsError(err)

    cutoff += 0.05

    print "DONE."
    
    # BLAST the nr with the seed sequences and download sequences
    seed_fasta = [s.formatFasta() for s in seed_list]
    localBlast(seed_fasta,"nr_blast.xml","nr/nr")
    seed_blast = parseBlastXML("nr_blast.xml",
                               tag_list=["Hit_def","Hit_id"])

    raw_homologs = [Homolog(s[0],s[1].split("|")[1].strip())
                    for s in seed_blast]
    cleaned_homologs = cleanHomologs(raw_homologs)

    to_download = [h.accession for h in cleaned_homologs]
    downloadSequences(to_download,"nr_download.xml")
    
    starting_homologs = parseFastaXML("nr_download.xml")
    starting_homologs = cleanHomologs(starting_homologs,min_length=min_length,
                                      max_length=max_length)
    
    # Make a low-redundancy set of sequences to reverse-blast starting 
    # homologs against.  Since the seeds are in this dataset, they will be
    # preferentially kept over any blast results because of their low ranks.
    # Since we are using an empirically-determined cutoff that will never
    # combine seed sequences, this will maximize the loss of redundancy, while
    # making sure we don't lose our seeds.  
    subject_homologs = runCdhit(starting_homologs+seed_list,redund_cutoff=cutoff) 
 
    # Reverse blast homologs within the species against the seed sequences.  
    seq2blastdb([h.formatFasta() for h in subject_homologs],
                db_name="reverse_query")

    f = open("coupling-list.txt","w")
    to_reverse_blast = seed_list + starting_homologs
    homolog_defs = dict([(h.unique_name,h.definition)
                         for h in to_reverse_blast])
    homolog_ranks = dict([(h.unique_name,h.rank)
                         for h in to_reverse_blast])

    counter = 0
    for h in to_reverse_blast:

        print h.organism, h.unique_name, h.definition
        localBlast([h.formatFasta()],"oB-tmp_result.xml","reverse_query",
                   force=True,hitlist_size=(num_homolog_cutoff+1))
        reverse_homologs = parseBlastXML("oB-tmp_result.xml")

        # Maybe we don't even find an ortholog?
        if len(reverse_homologs) == 0:
            top_hit = "NA"
            continue
    

        # Record the accession number and e value of the top num_homolog_cutoff
        # hits.

        f.write("%i\t%s\t%s" % (counter,h.definition,h.rank))
        
        num_homolog = 0
        index = 0
        while num_homolog < num_homolog_cutoff:
            if reverse_homologs[index][0] != h.unique_name:
                f.write("\t%s\t%s\t%s" % (homolog_defs[reverse_homologs[index][0]],
                                          homolog_ranks[reverse_homologs[index][0]],
                                          reverse_homologs[index][-1]))
                num_homolog += 1
            index += 1

        f.write("\n")
        f.flush()
        counter += 1
    
    f.close()
    

def main(argv=None):
    """
    Run ortholog blast!
    """

    # Parse the command line
    if argv == None:
        argv = sys.argv[1:]

    try:
        seed_fasta = argv[0]
    except IndexError:
        err = "Incorrect command line arguments!\n"
        err += "USAGE:\n\n%s\n\n" % __usage__
        raise BlastToolsError(err)


    reverseBlast(seed_fasta)



# Run main if we're starting from the command line
if __name__ == "__main__":
    main()        
