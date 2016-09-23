from __future__ import absolute_import
from phylogenetics.utils import SubclassError

def file(function):
    """Decorator to read data from a file.
    """
    def wrapper(self, data=None, path=None, *args, **kwargs):
        """ """
        # If a filename is not given, return string.
        if fname is None:
            # Make sure data was given since no file is named.
            if data is None:
                raise Exception("""`data` cannot be type, NoneType.""")
        # Write to a file
        else:
            try:
                # Try to write a straight string.
                with open(fname, "r") as f:
                    data = f.read()
            except:
                # IF writing string failed, try writing bytes.
                with open(fname, "rb") as f:
                    data = f.read()
        # Create a string of whatever datatype
        string = function(self, data, *args, **kwargs)
        return string
    return wrapper
