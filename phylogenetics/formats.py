# Easily convert between different phylogenetics file formats.
#
# Authors: Dr. Mike Harms
#          Zach Sailer
# -------------------------------------------------
# Fasta conversions
# -------------------------------------------------

import re

from phylogenetics.base import Homolog, HomologSet

# XML parser import
from xml.etree import ElementTree as ET

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
    num_columns = list(lengths.keys())[0]

    to_write = ["%s\n%s\n" % (o[:10],o[10:]) for o in split_out]
    to_write = "".join(to_write)

    final_out = "%i  %i\n\n%s\n" % (num_seq,num_columns-10,to_write)

    return final_out


# ---------------------------------------------------
# BLAST XML/fasta format
# ---------------------------------------------------

def flatten_concatenated_XML(input_file,key_tag):
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


def parse_blast_XML(filename,tag_list=("Hit_def","Hit_id")):
    """
        Parse XML file of hits returned after BLASTing a sequence
        against a database.

        Args:
        ----
        filename: str
            XML filename returned from BLAST
        tag_list: tuple
            Tuple of XML tags that will be included in output data for
            each sequence.

        Returns:
        -------
        all_hits: list of dicts
            List of sequence data for all hits in XML file (with key,values
            given by tag_list).
    """

    # Fix screwed up XML if blasts were done in series...
    blast_input = flatten_concatenated_XML(filename,"BlastOutput_iterations")

    # Read blast file properties (in tag_list) into a list to dump out
    blast_input = ET.XML(blast_input)

    all_hits = []

    # Navigate to proper level in XML
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

                    # Separate the tag_list out from hits and hit specs.
                    Hit_list = [t[4:] for t in tag_list if t[0:3] == "Hit"]
                    Hsp_list = [t[4:] for t in tag_list if t[0:3] == "Hsp"]

                    # Construct a dictionary of tags stripped from XML
                    properties = dict([(h.tag[4:],str(h.text).strip()) for h in hit])
                    subset = dict((k, properties[k]) for k in Hit_list)
                    all_hits.append(subset)

                    for property in hit:

                        if property.tag == "Hit_hsps":

                            for hsp in property:

                                # Construct a dictionary of tags stripped from XML
                                hsp_properties = dict([(p.tag[4:],str(p.text).strip())
                                                       for p in hsp])
                                subset = dict((k, hsp_properties[k]) for k in Hsp_list)

                                # Append inner Hsp properties to list
                                all_hits[-1].update(subset)
                                break

    return all_hits


def parse_blast_fasta(filename):
    """
        Parse Blast's Fasta/XML formatted file, returning a list of each
        sequence's data.

        Args:
        ----------
        filename: str
            Fasta/XML formatted file from Blast Output.

        Returns:
        -------
        sequences: list of dicts
            List of sequences data in their own lists.
    """

    # Fix screwed up XML because sequences downloaded and output concatenated
    sequence_input = flatten_concatenated_XML(filename,"TSeqSet")

    # Now we should have valid XML...
    sequence_input = ET.XML(sequence_input)

    # XML Tag prefix to strip
    prefix = "TSeq_"

    sequences = []
    for i, sequence in enumerate(sequence_input):
        # Rip out all properties of a sequence
        properties = dict([(s.tag[len(prefix):],str(s.text).strip()) for s in sequence])

        # Append to sequences.
        sequences.append(properties)

    return sequences


def blast_to_homologset(filename, tag_list=()):
    """ Load blast XML file as HomologSet. """

    # Don't discriminate on the type of Blast XML file format, try both.
    try:
        # First, try blasting hits format.
        hits = parse_blast_XML(filename, tag_list=tag_list)
    except:
        # Then, try blast fasta format
        hits = parse_blast_fasta(filename)

    homologs = []
    for i in range(len(hits)):
        # Make a unique id for each sequence.
        unique_id = "XX%08d" % i

        # Create homolog instance for each sequence.
        homologs.append(Homolog(unique_id, **hits[i]))

    # Build homolog set from xml file
    homologset = HomologSet(homologs)
    return homologset
