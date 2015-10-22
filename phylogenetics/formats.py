# Easily convert between different phylogenetics file formats.
#
# Authors: Dr. Mike Harms
#          Zach Sailer
# -------------------------------------------------
# Fasta conversions
# -------------------------------------------------


# ----------------------------------------------------
# Methods for parsing different file formats
# ----------------------------------------------------

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
