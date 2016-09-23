def file(function):
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
