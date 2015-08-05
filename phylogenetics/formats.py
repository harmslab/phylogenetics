# Easily convert between different phylogenetics file formats.s
#
# -------------------------------------------------
# Fasta conversions
# -------------------------------------------------

from phylogenetics.base import Homolog, HomologSet

class Fasta2PhylipError(Exception):
    """
    General error class for this module.
    """

    pass


def fasta2phylip(lines):
    """
    Take the lines from a fasta file and return in fasta format.
    """

    lines = [l for l in lines if l.strip() != ""]

    # Read fasta file
    out = []
    for l in lines:
        if l[0] == ">":
            if len(out) != 0:
                out.append("\n")
            out.append("%-10s" % l[1:11])
        else:
            out.append(l.strip())

    out = "".join(out)
    split_out = out.split("\n")

    # Quick sanity check
    lengths = dict([(len(l),i) for (i,l) in enumerate(split_out)])
    if len(lengths) > 1:
        err = "Some sequences have different lengths!\n"
        for d in lengths.keys():
            err = "%s, %i\n" % (split_out[lengths[d]][:10],d)

        raise Fasta2PhylipError(err)

    num_seq = len(split_out)
    num_columns = lengths.keys()[0]

    to_write = ["%s\n%s\n" % (o[:10],o[10:]) for o in split_out]
    to_write = "".join(to_write)

    final_out = "%i  %i\n\n%s\n" % (num_seq,num_columns-10,to_write)

    return final_out
    
    
# ---------------------------------------------------
# BLAST XML/fasta format
# ---------------------------------------------------
    
def load_blast_xml(filename):
    """ Load blast XML file as homolog objects. """

    def substring(s, first, last):
        """ Function for returning a substring between two tags
        
            Example:
            -------
            s = "letmeshowyousomethingcool"
            first = "letme"
            last = "somethingcool"
            
            Returns 
            >>> "showyou", 0, 11
            
        """
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end], start, end

    # load file
    f = open(filename, "r")
    string = f.read()
    f.close()
    
    jump = 0
    homologs = list()
    
    while jump < len(string):
        # Search for the next instance of <TSeq> tag for sequence data
        try:
            # Find individual sequence data
            sequence_data, start, end = substring(string[jump:], "<TSeq>", "</TSeq>")
        
            # Find all sequence info tags
            seq_tags = list()
            sub_jump = 0
        
            # Find all unique tag in sequence data.
            while jump < len(sequence_data):
                try:
                    tag, sub_start, sub_end = substring(sequence_data[sub_jump:], "<TSeq_", ">")
                    sub_jump += sub_end
                    seq_tags.append(tag)
                # Once no one new tags are found, break loop
                except:
                    break
        
            # Find and strip sequence data from sequence tags found previously
            kwargs = {}
            # ignore first xml tag, TSeq_seqtype
            for t in seq_tags[1:]:
                # Get data, x and y are junk for this purpose
                data, x,y = substring(sequence_data, "<TSeq_" + t + ">", "</TSeq_" + t + ">")
                kwargs[t] = data
        
            # Make a unique ID for this sequence
            unique_id = "XX%08d" % counter
        
            # Add this homolog to a set of homologs
            h = Homolog(unique_id, **kwargs)
            homologs.append(h)
        
            # Update jump
            jump += end
        # Once no more tags are found, break loop
        except:
            break
            
    return homologs