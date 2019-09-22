from time import localtime, strftime
import os
import json
from functools import wraps

def track_in_history(method):
    """Track initialization, changes to data, and calculations in a history
       json."""
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        """"""
        # Run the method
        output = method(self, *args, **kwargs)

        # Get the time
        time =strftime("%Y-%m-%d %H:%M:%S", localtime())

        # Create a history item
        history = {'time':time, 'method':str(method), 'args':str(args), 'kwargs':kwargs}

        # Append to main history list
        self.history.append(history)

        # Write out current dataframe to csv
        #new_df.to_csv(self.project_dir + "/" + self.history[0]["method"].split(".")[1].split(" ")[0]) + ".csv"

        # Write history to a json file.
        history_file = os.path.join(self.project_dir, 'history.json')
        with open(history_file, 'w') as f:
            json.dump(self.history, f)

        return output
    return wrapper
