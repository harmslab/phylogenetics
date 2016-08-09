from .utils import get_time

def BaseClass(object):
    """Provides a template class for adding/removing attributes and metadata.
    """
    def __init__(self):
        self._metadata = {}
        # On initialization, set the datetime of the
        self.addattr(datetime=get_time())

    @property
    def metadata(self):
        """Get metadata for this object"""
        return self._metadata

    def addattr(self, **kwargs):
        """Add attributes to object and metadata."""
        for key, value in kwargs:
            setattr(self, key, value)
            self._metadata[key] = value

    def rmattr(self, *arg):
        """Remove attribute from object and metadata."""
        for a in args:
            delattr(self, a)
            del self._metadata[a]
