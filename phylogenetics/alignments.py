import warnings

from . import handlers
from .sequences import Sequence, SequenceList

class AlignedSequence(handlers.Handler):
    """Aligned sequence.
    """
    def __init__(self, sequence, **kwargs):
        super(AlignedSequence, self).__init__(sequence=sequence, **kwargs)

    @property
    def _prefix(self):
        return "Ali"

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
        self.addattr(Sequence=id)
        self.Sequence = Sequence

    @handlers.history
    def unlink(self):
        """Unlink a sequence to this Alignment.
        """
        # Try to remove a Homolog
        try:
            self.rmattr(Sequence)
            delattr(self, "Sequence")
        except AttributeError:
            warnings.warn("No Sequence was linked to the AlignedSequence object.")


class Alignment(handlers.HandlerContainer):
    """Alignment object. Contains a collection of AlignedSequence.
    """
    def __init__(self, *AlignedSequences, **kwargs):
        super(Alignment, self).__init__(*AlignedSequences, **kwargs)

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
        for Sequence in SequenceList.list:
            obj = type(Sequence)
            if obj in self._link_types:
                id_number = Sequence.id[3:]
                aligned_seq = self._contents["Ali"+id_number]
                aligned_seq.link(Sequence)
            else:
                raise Exception

    @handlers.history
    def unlink(self, *ids):
        """Link Sequence from Alignment
        """


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
