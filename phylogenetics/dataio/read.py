# Module for reading HomologSets
import json
import pickle

def read_from_file(function):
    """ Decorator to read data from a file.
    """
    def wrapper(self, fname=None, data=None, *args, **kwargs):
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
                with open(fname, "w") as f:
                    data = f.read(string)
            except:
                # IF writing string failed, try writing bytes.
                with open(fname, "wb") as f:
                    data = f.read(string)
            return string

        # Create a string of whatever datatype
        string = function(self, data, *args, **kwargs)

    return wrapper

class Read(object):

    def __init__(self, HomologSet):

        self._homologset = HomologSet

    @read_from_file
    def alignment(self, data):
        """ """
        pass

    @read_from_file
    def fasta(self, data):
        """ Read a fasta string.

        """
        pass


    @read_from_file
    def csv(self, data):
        pass

    @read_from_file
    def json(self, data):
        pass

    @read_from_file
    def phylip(self, data):
        pass

    @read_from_file
    def newick(self, data):
        pass
