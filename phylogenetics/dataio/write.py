from __future__ import absolute_import

def file(function):
    """ Decorator to write the output of a function to a file if `fname` kwarg exists.
    Otherwise, output is just returned.

    Output must be a string.
    """
    def wrapper(self, path=None, *args, **kwargs):
        """ """
        # Create a string of whatever datatype
        string = function(self, *args, **kwargs)
        # If a filename is not given, return string.
        if path is None:
            return string
        # Write to a file
        else:
            try:
                # Try to write a straight string.
                with open(path, "w") as f:
                    f.write(string)
            except:
                # IF writing string failed, try writing bytes.
                with open(path, "wb") as f:
                    f.write(string)
            #return string
    return wrapper
