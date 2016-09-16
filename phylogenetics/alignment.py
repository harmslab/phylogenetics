from . import handlers

class AlignedHomolog(handlers.Handler):
    """Aligned sequence.
    """
    def __init__(self, **kwargs):
        super(AlignedHomolog, self).__init__(**kwargs)

    def _link_Homolog(self, Homolog):
        """Link a Homolog object to this Alignment
        """
        pass

class Alignment(handlers.ContainerHandler):
    """Alignment object. Contains a collection of AlignedHomologs.
    """
    def __init__(self, *AlignedHomologs, **kwargs):
        super(Alignment, self).__init__(*AlignmentHomologs, **kwarg)

    @property
    def alignment(self):
        return self._contents

    @property
    def _prefix(self):
        return "ALI"

    @property
    def _child_type(self):
        return AlignedHomolog

    def _link_Homologs(self, *Homologs):
        """Link Homologs to Alignment
        """
        pass
