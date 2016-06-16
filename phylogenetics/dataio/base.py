"""Base module for writing
"""
from phylogenetics.utils import SubclassError

class Read(object):

    def _data_to_object(self):
        """Read object from dictionary from pickled, json-ed, etc."""
        raise SubclassError("""Must be implemented in subclass""")

    def _sequences_to_object(self):
        """Read object from list of tuples, from formats like fasta,
        phylip, etc.
        """
        raise SubclassError("""Must be implemented in subclass""")

class Write(object):

    def _object_to_data(self):
        """Write object to a dictionary for pickling, json, etc."""
        raise SubclassError("""Must be implemented in subclass""")

    def _object_to_sequences(self):
        """Write object to list of tuples for formats such as fasta,
        phylip, etc.
        """
        raise SubclassError("""Must be implemented in subclass""")


def read_from_file(function):
    """ Decorator to read data from a file.
    """
    def wrapper(self, data=None, fname=None, *args, **kwargs):
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


def write_to_file(function):
    """ Decorator to write the output of a function to a file if `fname` kwarg exists.
        Otherwise, output is just returned.

        Output must be a string.
    """
    def wrapper(self, fname=None, *args, **kwargs):
        """ """
        # Create a string of whatever datatype
        string = function(self, *args, **kwargs)

        # If a filename is not given, return string.
        if fname is None:
            return string

        # Write to a file
        else:
            try:
                # Try to write a straight string.
                with open(fname, "w") as f:
                    f.write(string)
            except:
                # IF writing string failed, try writing bytes.
                with open(fname, "wb") as f:
                    f.write(string)
            #return string

    return wrapper
