from __future__ import absolute_import

import json

def read(data):
    """
        Read pickle string and convert to object.
    """
    object_data = json.loads(data)
    return object_data

def write(object_data):
    """
        Write data as pickle string.
    """
    output = json.dumps(object_data)
    return output
