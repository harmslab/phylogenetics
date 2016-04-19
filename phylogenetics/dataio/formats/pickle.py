import pickle

def write(object_data):
    """
        Write data as pickle string.
    """
    output = pickle.dumps(object_data)
    return output


def read(data):
    """
        Read pickle string and convert to object.
    """
    object_data = pickle.loads(data)
    return object_data
