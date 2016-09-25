from . import handlers

class Sequence(handlers.Handler):
    """Creates an object for managing data for a single homolog sequence.

    Parameters
    ----------
    sequence : string
        Sequence as a string.

    Keyword arguments
    -----------------
    """
    def __init__(self, sequence, **kwargs):
        super(Sequence, self).__init__(sequence=sequence, **kwargs)

    @property
    def _prefix(self):
        return "Seq"

class SequenceList(handlers.HandlerContainer):
    """A Sequence container object.
    """
    def __init__(self, *Sequences, **kwargs):
        super(SequenceList, self).__init__(*Sequences, **kwargs)

    @property
    def sequences(self):
        """Return dictionary with sequence ids mapped to Sequence objects."""
        return self._contents

    @property
    def _prefix(self):
        return "Sequencelist"

    @property
    def _child_types(self):
        return [Sequence]
