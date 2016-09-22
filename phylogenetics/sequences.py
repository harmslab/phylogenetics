from . import handlers
from .dataio import sequenceio, sequencelistio

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
        # Attach Write-ing object
        #self.Read = homologio.Read(self)
        self.write = sequenceio.Write(self)
        self.read = sequenceio.Read(self)

    @property
    def _prefix(self):
        return "Seq"

class SequenceList(handlers.HandlerContainer):
    """A homolog container object.
    """
    def __init__(self, *Sequences, **kwargs):
        super(SequenceList, self).__init__(*Sequences, **kwargs)
        self.write = sequencelistio.Write(self)
        self.read = sequencelistio.Read(self)
        #self.Write = sequencelistio.Write(self)

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
