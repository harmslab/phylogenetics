# Easily convert between different phylogenetics file formats.
#
# Authors: Dr. Mike Harms
#          Zach Sailer
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

def flattenConcatenatedXML(input_file,key_tag):
    """
    Clean up naively concatenated XML files by deleting begin/end tags that
    occur at the place where the two files were concatenated.
    NOTE: This will break and break royally if the key_tags are on the same
    lines as other important entries.
    """

    f = open(input_file,'r')
    input = f.readlines()
    f.close()

    set_start = re.compile("<%s>" % key_tag)
    set_end =   re.compile("</%s>" % key_tag)

    # Find all beginning and end tags...
    starts = [i for i, l in enumerate(input) if set_start.search(l) != None]

    # If this tag occurs more than once...
    if (len(starts) != 1):

        # Keep the first start reverse so we are chewing from the bottom.
        starts.pop(0)
        starts.reverse()

        # Remove all lines between each end and start, again chewing from the
        # bottom.
        for i in range(len(starts)):
            e = starts[i]
            while set_end.search(input[e]) == None:
                input.pop(e),
                e = e - 1
            input.pop(e)

    # Return freshly minted, clean XML
    return "".join(input)

def parseBlastXML(blast_file,tag_list=("Hit_def","Hit_id")):
    """
    Parse BLAST xml output, extracting tags specified in tag_list and putting
    into a list.  E-value is always appended after the last requested tag.
    """

    # Fix screwed up XML if blasts were done in series...
    blast_input = flattenConcatenatedXML(blast_file,"BlastOutput_iterations")

    # Read blast file properties (in tag_list) into a list to dump out
    blast_input = ET.XML(blast_input)

    all_hits = []
    for blast_output in blast_input:
        if blast_output.tag != "BlastOutput_iterations":
            continue

        for iteration in blast_output:
            if iteration.tag != "Iteration":
                continue

            for hits in iteration:
                if hits.tag != "Iteration_hits":
                    continue

                for hit in hits:
                    Hit_list = [t for t in tag_list if t[0:3] == "Hit"]
                    Hsp_list = [t for t in tag_list if t[0:3] == "Hsp"]

                    properties = dict([(h.tag,str(h.text)) for h in hit])
                    all_hits.append([properties[t] for t in Hit_list])

                    for property in hit:
                        if property.tag == "Hit_hsps":
                            for hsp in property:
                                hsp_properties = dict([(p.tag,str(p.text))
                                                       for p in hsp])

                                # Append inner Hsp properties to list
                                all_hits[-1] += ([hsp_properties[t] for t in Hsp_list])
                                break

    return all_hits

    
def Blast_to_homologs(filename):
    """ Load blast XML file as homolog objects. """







def Blast_to_homologs(filename):
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

    jump, counter = 0, 0
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
            while sub_jump < len(sequence_data):
                try:
                    tag, sub_start, sub_end = substring(sequence_data[sub_jump:], "<TSeq_", ">")
                    sub_jump += sub_end
                    seq_tags.append(tag)
                # Once no one new tags are found, break loop
                except ValueError:
                    sub_jump = len(sequence_data)

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
            counter += 1
        # Once no more tags are found, break loop
        except ValueError:
            jump = len(string)

    return homologs
