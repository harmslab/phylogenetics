__doc__ = """Defines handlers for all objects in the phylogenetics package.

Basic Handler template has a set of attributes that are written to disk.

Basic ContainerHandler template has contents with multple objects.
"""

import time, copy

def history(method):
    """Update the history of the object by adding the current time to history
    list.
    """
    def update_history(*args, **kwargs):
        now = time.strftime("%c")
        action = method.__name__
        obj = args[0]
        info = (now, action)
        obj.history.append(info)
        return method(*args, **kwargs)
    return method

def Handler(object):
    """Provides a template class for adding/removing attributes and metadata.
    """
    def __init__(self, **kwargs):
        self.attrs = {}
        # On initialization, set the datetime
        now = time.strftime("%c")
        action = "init"
        info = (now, action)
        self.history = [info]
        self.addattr(**kwargs)

    @history
    def addattr(self, **kwargs):
        """Add attributes to object and metadata."""
        for key, value in kwargs:
            setattr(self, key, value)
            self.attrs[key] = value

    @history
    def rmattr(self, *arg):
        """Remove attribute from object and metadata."""
        for a in args:
            delattr(self, a)
            del self.attrs[a]


def HandlerContainer(Handler):
    """A general container to manage a collection instances of a the same object.

    All changes to the object are stored in history.
    """
    def __init__(self, *arg, **kwargs):
        super(Container, self).__init__(**kwargs)
        self._contents = {}
        self.add(*args)

    @property
    def metadata(self):
        """Get the attributes and content metadata of this object"""
        data = dict([(c.id, c.metadata), for c in self._contents])
        data.update(**self.attrs)
        return data

    @property
    def _type_exception_message(self):
        return """And object already exists in contents with the same ID."""

    @property
    def _prefix(self):
        raise Exception("""Must be implemented in a Subclass.""")

    @property
    def _child_type(self):
        """The specific objects that are contained in this object.
        Note that child type must be a list. If only a single object,
        """
        raise Exception("""Must be implemented in a Subclass.""")

    def _check_type(self, item):
        """ Check that the item is an expected object.
        """
        if item.__class__ not in self._child_type:
            raise Exception("Argument must be a(n) `" + \
                self._child_type.__name__ + "` object!")

    def _assign_id(self, item):
        """Assigns an `id` to object, with a given prefix.
        """
        # Check that argument is an expected object
        self._check_type(item)
        # If the object doesn't already have an ID, give it one
        if hasattr(item, "id") is False:
            number = 0
            prefix = item._prefix
            num_size = 10 - len(prefix)
            label = str(number)
            new_id = prefix + "0"*(numsize - len(label)) + label
            while new_id in self._contents.keys():
                number += 1
                label = str(number)
                new_id = prefix + "0"*(numsize - len(label)) + label
            item.addattr(id=new_id)
        # If ID exists in object, check that it isn't in this object already
        else:
            if item.id in self._contents:
                raise Exception(self._type_exception_message)

    @history
    def add(self, *args):
        """Add an instance of the main object type to this container.
        """
        for a in args:
            self._assign_id(a)
            setattr(self, a.id, a)
            self._contents[a.id] = a

    @history
    def rm(self, *ids):
        """Remove an instance of the main object from this container.
        """
        for id in ids:
            delattr(self, id)
            del self._contents[id]
