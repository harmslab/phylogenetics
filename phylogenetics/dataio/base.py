# Base submodule for dataio module

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
