from . import handlers

class Ancestor(handlers.Handler):
    """Ancestor message
    """
    def __init__(self, **kwargs):
        super(Ancestor, self).__init__(**kwargs)

    @property
    def _prefix(self):
        return "Anc"

class AncestorList(handlers.HandlerContainer):
    """
    """
    def __init__(self, *Ancestors, **kwargs):
        super(AncestorList, self).__init__(*Ancestors, **kwargs)

    @property
    def _prefix(self):
        return "AncList"

    @property
    def _child_type(self):
        """Object to ancestor."""
        return [Ancestor]
