# XML parser import
from __future__ import absolute_import

import re
from xml.etree import ElementTree as ET

def read(xml_string):
    """Parse Blast's Fasta/XML formatted file, returning a list of each
    sequence's data.

    Parameters
    ----------
    xml_string : str
        Fasta/XML formatted string from Blast Output.

    Returns
    -------
    sequences : list of dicts
        List of sequences data in their own lists.
    """
    # Fix screwed up XML because sequences downloaded and output concatenated
    sequence_input = flatten_concatenated_XML(xml_string, "TSeqSet")

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


def flatten_concatenated_XML(xml_string,key_tag):
    """Clean up naively concatenated XML files by deleting begin/end tags that
    occur at the place where the two files were concatenated.
    NOTE: This will break and break royally if the key_tags are on the same
    lines as other important entries.
    """
    input = xml_string.split("\n")
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
