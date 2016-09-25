import warnings

from . import handlers
from .sequences import Sequence, SequenceList
from .utils import LinkError

class AlignedSequence(handlers.Handler):
    """Aligned sequence.
    """
    def __init__(self, sequence, **kwargs):
        super(AlignedSequence, self).__init__(links=[], sequence=sequence, **kwargs)

    @property
    def _prefix(self):
        return "Seq"

    @property
    def _link_types(self):
        """"""
        return [Sequence]

    @handlers.history
    def link(self, Sequence):
        """Link a Sequence object to this Alignment.

        A call-by-reference link is made to the Sequence, and a new attribute is
        added to self.
        """
        id = Sequence.id
        # Get old links
        links = self.links
        # append new link
        links.append(id)
        # update link attribute
        self._addattr(links=links)
        # add Sequence to object attribute, without included in contents.
        self.Sequence = Sequence

    @handlers.history
    def unlink(self):
        """Unlink a sequence to this Alignment.
        """
        # Try to remove a Homolog
        try:
            delattr(self, "Sequence")
            self._addattr(links=[])
        except AttributeError:
            warnings.warn("No Sequence was linked to the AlignedSequence object.")


class Alignment(handlers.HandlerContainer):
    """Alignment object. Contains a collection of AlignedSequence.
    """
    def __init__(self, *AlignedSequences, **kwargs):
        super(Alignment, self).__init__(links=[], *AlignedSequences, **kwargs)

    @property
    def alignment(self):
        return self._contents

    @property
    def _prefix(self):
        return "Alignment"

    @property
    def _child_types(self):
        """"""
        return [AlignedSequence]

    @property
    def _link_types(self):
        """"""
        return [SequenceList]

    @handlers.history
    def link(self, SequenceList):
        """Link SequenceList to Alignment.
        """
        # Check that a SequenceList was given
        if type(SequenceList) in self._link_types:
            # Get links attribute
            links = self.links
            # Update attribute
            links.append(SequenceList.id)
            self.addattr(links=links)
            # Add Sequencelist as attribute not stored in metadata.
            self.SequenceList = SequenceList
            # Link individual sequences too.
            for AlignedSequence in self.list:
                id = AlignedSequence.id
                Sequence = self.SequenceList.contents[id]
                AlignedSequence.link(Sequence)
        else:
            raise LinkError("Argument must be type == SequenceList.")

    @handlers.history
    def unlink(self):
        """Remove SequenceList from Alignment
        """
        # Try to remove a SequenceList
        try:
            # Remove links from individual aligned seqs.
            for aligned_seq in self.list:
                aligned_seq.unlink()
            # Remove sequencelist if it exists.
            delattr(self, "SequenceList")
            self._addattr(links=[])
        except AttributeError:
            warnings.warn("No Sequence was linked to the AlignedSequence object.")

class AlignmentList(handlers.HandlerContainer):
    """Container object for managing multiple Alignments.
    """
    def __init__(self, *Alignments, **kwargs):
        super(AlignmentList, self).__init__(*Alignments, **kwargs)

    @property
    def _prefix(self):
        return "AlignmentList"

    @property
    def alignments(self):
        return self._contents

    @property
    def _child_types(self):
        return [Alignment]
