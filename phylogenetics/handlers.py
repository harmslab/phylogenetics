__doc__ = """Basic ``Handler`` objects for managing phylogenetic data.
"""
import time, copy, json, importlib
from .dataio import read, write, formats

def history(method):
    """Update the history of the object by adding the current time to history
    list.
    """
    def update_history(*args, **kwargs):
        now = time.strftime("%c")
        action = method.__name__
        obj = args[0]
        info = [now, action]
        obj.history.append(info)
        obj.attrs["history"].append(info)
        return method(*args, **kwargs)
    return update_history

class Handler(object):
    """Provides a template class that attaches an API to metadata dictionary.

    Example
    -------
    ..code::

        metadata = {
            "type" : "MyObject",
            "module" : "foo.bar",
            "history" = ["Tues September 27 9:00:00AM", "init"]
        }

    
    """
    def __init__(self, **kwargs):
        # Some basic attributes to save
        self.attrs = {}
        self._addattr(
            type=type(self).__name__,
            module=self.__module__
        )
        # Initialize history keywords.
        self._init_history(**kwargs)
        self._addattr(**kwargs)

    def _init_history(self, **kwargs):
        """Initialize history."""
        if "history" not in kwargs:
            # On initialization, set the datetime
            now = time.strftime("%c")
            action = "init"
            info = [now, action]
            self.attrs["history"] = [info]
        else:
            self.attrs["history"] = kwargs["history"]
        self.history = self.attrs["history"]

    def print(self):
        """Print metadata as formated json."""
        print(json.dumps(self.metadata, indent=4, separators=(',', ': ')))

    @property
    def _schemas(self):
        """"""
        return Exception("""Must be implemented in a subclass.""")

    @property
    def _prefix(self):
        """"""
        raise Exception("""Must be implemented in a Subclass.""")

    @property
    def metadata(self):
        """alias for attrs"""
        return self.attrs

    def _metadata_to_object(self, metadata):
        self._addattr(**metadata)

    @history
    def addattr(self, **kwargs):
        """Add attributes to object and metadata."""
        self._addattr(**kwargs)

    def _addattr(self, **kwargs):
        """see addattr doc."""
        for key, value in kwargs.items():
            setattr(self, key, value)
            self.attrs[key] = value

    @history
    def rmattr(self, *arg):
        """Remove attribute from object and metadata."""
        self._rmattr(*arg)

    def _rmattr(self, *args):
        """see rmattr doc."""
        for a in args:
            delattr(self, a)
            del self.attrs[a]

    @history
    def append_to_attr(self,**kwargs):
        """Update an attribute of the Handler. Only work for arguments that are
        lists or dictionarys. If Argument is passed, append to list of
        """
        for key, value in kwargs.items():
            old = getattr(self, key)
            old.append(value)
            self.attrs[key] = old
            setattr(self, key, old)

class HandlerContainer(Handler):
    """A general container to manage a collection instances of a the same object.

    Note
    ----
    All changes to the object are stored in the history attribute.
    """
    def __init__(self, *args, **kwargs):
        super(HandlerContainer, self).__init__(**kwargs)
        self._contents = {}
        self._add(*args)

    @classmethod
    def get(cls, data, schema, **kwargs):
        """"""
        instance = cls()
        instance.read(data=data, schema=schema)
        return instance

    @write.file
    def write(self, schema):
        """Write Container to metadata."""
        writer = getattr(formats, schema).write
        try:
            return writer(self.metadata)
        except:
            contents = self.metadata["contents"]
            return writer(contents)

    @history
    @read.file
    def read(self, data, schema):
        """Read data with schema to object."""
        parser = getattr(formats, schema).read
        metadata = parser(data)
        self._metadata_to_object(metadata)

    #@read.file
    def update(self, obj, data, schema):
        """Update one of the pieces of a phylogenetics project which already
        contains data.
        """

    def _metadata_to_object(self, metadata):
        """parse metadata and construct objects from its info.
        """
        contents = metadata.pop("contents")
        # Add contents.
        for m in contents:
            mod = importlib.import_module(m["module"])
            Obj = getattr(mod, m["type"])
            # Initialize the object
            handler = Obj()
            # If its also a container, add its contents
            handler._metadata_to_object(m)
            self._add(handler)
        # Add other attributes to objects
        self._addattr(**metadata)

    @property
    def contents(self):
        return self._contents

    @property
    def list(self):
        """A list of the Handlers contained in this Container."""
        return [content for content in self._contents.values()]

    @property
    def metadata(self):
        """Get the attributes and content metadata of this object"""
        data = {"contents": [item.metadata for item in self.list]}
        data.update(**self.attrs)
        return data

    @property
    def _type_exception_message(self):
        return """And object already exists in contents with the same ID."""

    @property
    def _child_types(self):
        """The specific objects that are contained in this object.
        Note that child type must be a list. If only a single object,
        """
        raise Exception("""Must be implemented in a Subclass.""")

    def _check_type(self, item):
        """ Check that the item is an expected object.
        """
        if item.__class__ not in self._child_types:
            names = [c.__name__ for c in self._child_types]
            raise Exception("Argument must be one of the following objects: " + ", ".join(names))

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
            new_id = prefix + "0"*(num_size - len(label)) + label
            while new_id in self._contents.keys():
                number += 1
                label = str(number)
                new_id = prefix + "0"*(num_size - len(label)) + label
            item._addattr(id=new_id)
        # If ID exists in object, check that it isn't in this object already
        else:
            if item.id in self._contents:
                raise Exception(self._type_exception_message)

    @history
    def add(self, *args):
        """Add an instance of the main object type to this container.
        """
        self._add(*args)

    def _add(self, *args):
        """see add doc.
        """
        for a in args:
            self._assign_id(a)
            setattr(self, a.id, a)
            self._contents[a.id] = a

    @history
    def rm(self, *ids):
        """Remove an instance of the main object from this container.
        """
        self._rm(*ids)

    def _rm(self, *ids):
        """See rm doc.
        """
        for id in ids:
            delattr(self, id)
            del self.content[id]
