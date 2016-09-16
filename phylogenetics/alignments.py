import warnings

from . import handlers
from .sequences import Sequence

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
        """Link a sequence object to this Alignment
        """
        id = Sequence.id
        self.addattr(homolog=id)
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
            warnings.warn("No Homolog was linked to the AlignedSequence object.")


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
        return [AlignedSequence]

    @handlers.history
    def link(self, *Sequences):
        """Link Sequences to Alignment
        """
        for Sequence in Sequences:
            id_number = Sequence.id[3:]
            aligned_seq = self._contents["Ali"+id_number]
            aligned_seq._link_sequence(Sequence)

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
