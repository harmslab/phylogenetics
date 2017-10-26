import pandas as pd
from time import localtime, strftime

def track_in_history(method):
    """Track this call in the history DataFrame"""
    
    def wrapper(self, *args, **kwargs):
        
        # Prepare data for history dataframe
        time =strftime("%Y-%m-%d %H:%M:%S", localtime())
        args_as_str = ",".join([str(a) for a in args])
        kwargs_as_str = ",".join([str((key, val)) for key, val in kwargs.items()])
        
        # Create row
        history = pd.DataFrame({'time':[time], 
            'method':[method.__name__], 
            'args':[args_as_str], 
            'kwargs':[kwargs_as_str]}, dtype=str)

        # Append to main history dataframe
        self.history = self.history.append(history, ignore_index=True)
        
        # Now run method
        return method(self, *args, **kwargs)
    return wrapper

class Project(object):
    """
    """
    def __init__(self, project_dir):
        time = strftime("%Y-%m-%d %H:%M:%S", localtime())
        self.project_dir = project_dir
        self.history = pd.DataFrame({'time':[time], 
            'method':['__init__'], 
            'args':[project_dir], 
            'kwargs':[None]}, dtype=str)
         
    @track_in_history
    def test(self, blah, blah2=None):
        return 'hellos'
