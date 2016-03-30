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

class HomologFormats:

    def __init__(self, homolog):
        self.homolog = homolog

    def fasta(self, tags=None, aligned=False):
        """ Return fasta formatted string with named tags (in order given).

            If no tags are given, prints id with sequence.
        """
        f = ">"
        if tags is not None:
            # Remove sequence if its an elements
            if "sequence" in tags:
                tags.remove("sequence")
            f += "|".join([str(getattr(self.homolog, t)) for t in tags])
        else:
            f += self.homolog.id

        # New line with full sequence string (aligned if told to)
        if aligned:
            f += "\n" + self.homolog.latest_align +"\n"
        else:
            f += "\n" + self.homolog.sequence +"\n"

        return f

    def phylip(self, tags=None, **kwargs):
        """ Return a PhyML formatted string to write to file.

            Only allowed when Homolog has alignment attribute
        """
        # Check that an alignment attribute exists in homolog.
        if hasattr(self.homolog, "latest_align") is False:
            raise Exception("Homolog must have an alignment attribute to write phylip.")

        # Use specified alignment key
        if tags is not None:
            alignment = getattr(self.homolog, tags)
        else:
            alignment = self.homolog.latest_align

        f = "%s\n%s\n" % (self.homolog.id, alignment)
        return f

    def json(self, **kwargs):
        """ Return json formatted string. """
        return json.dumps(self.homolog.__dict__)

    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self.homolog)

    def csv(self, tags=None, header=True, delimiter=",", **kwargs):
        """ write csv string. """
        # Get all attributes if tags are not specified
        if tags is None:
            tags = list(vars(self.homolog).keys())

        # Add header is specified
        if header:
            f = delimiter.join(tags)
        else:
            f = ""

        # Try to print tag. If it doesnt exist, then leave blank
        vals = list()
        for t in tags:
            try:
                vals.append(str(getattr(self.homolog,t)))
            except:
                vals.append("")
        f += delimiter.join(vals) + "\n"

        return f

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
