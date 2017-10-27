import pandas
import toytree
import phylopandas
import dendropy

from time import localtime, strftime

def track_in_history(method):
    """Track this call in the history DataFrame"""
    def wrapper(self, *args, **kwargs):
        """"""
        # Prepare data for history dataframe
        time =strftime("%Y-%m-%d %H:%M:%S", localtime())
        args_as_str = ",".join([str(a) for a in args])
        kwargs_as_str = ",".join([str((key, val)) for key, val in kwargs.items()])
        
        # Create row
        history = pandas.DataFrame({'time':[time], 
            'method':[method.__name__], 
            'args':[args_as_str], 
            'kwargs':[kwargs_as_str]}, dtype=str)

        # Append to main history dataframe
        self.history = self.history.append(history, ignore_index=True)
        
        # Now run method
        return method(self, *args, **kwargs)
    return wrapper

class TreeProject(object):
    """Main Phylogenetics Project class."""
    def __init__(self, project_dir):
        time = strftime("%Y-%m-%d %H:%M:%S", localtime())
        self.data = {'leafs':None, 'ancestors':None, 'tree':None}
        self.project_dir = project_dir
        self.history = pandas.DataFrame({'time':[time], 
            'method':['__init__'],
            'args':[project_dir], 
            'kwargs':[None]}, dtype=str)
          
        
    def __str__(self):
        # Get history
        history = self.history.iloc[-1]
        leafs_bool = self.data['leafs'] is not None
        ancestors_bool = self.data['ancestors'] is not None
        tree_bool = self.data['tree'] is not None
        
        # Start building history.
        info = ["TreeProject(project_dir={})\n".format(self.project_dir),
                "    last modified\t{}\n".format(history['time']),
                "    last edit\t\t{}\n".format(history['method'])]
                
        # If leafs, add stats
        info += ["    leafs\t\t{}\n".format(True)]
        if leafs_bool:
            leafs = self.data['leafs']
            info += ["        - num of leafs\t{}\n".format(len(leafs))]
    
        # If ancestors, add stats.
        info += ["    ancestors\t\t{}\n".format(True)]
        if ancestors_bool:            
            ancs = self.data['leafs']
            info += ["        - num of ancs\t{}\n".format(len(ancs))]            
        
        # If tree, add stats.
        info += ["    tree\t\t{}\n".format(True)]
        return "".join(info).strip()
                
    def __repr__(self):
        return self.__str__()
         
    @track_in_history
    def test(self, blah, blah2=None):
        return self

    def _add_leafs(self, data):
        """"""
        if isinstance(data, pandas.DataFrame) == False and isinstance(data, phylopandas.DataFrame) == False:
            raise Exception('Bad datatype.')
        self.data['leafs'] = data
        self.leafs = data

    def _add_ancestors(self, data):
        """"""
        if isinstance(data, pandas.DataFrame) == False and isinstance(data, phylopandas.DataFrame) == False:
            raise Exception('Bad datatype.')
        self.data['ancestors'] = data
        self.ancestors = data
        
    def _add_tree(self, data):
        """"""
        if isinstance(data, dendropy.Tree) == False:
            raise Exception('Bad datatype.')
        self.data['tree'] = data
        self.tree = data

    @track_in_history
    def add_data(self, dtype, data):
        """
        """
        if dtype not in ['leafs', 'ancestors', 'tree']:
            raise Exception('dtype is not valid.')
        
        # Call method on data.
        method = getattr(self, '_add_{}'.format(dtype))
        method(data)
        return self

    @track_in_history
    def read_data(self, dtype, path, schema='fasta', **kwargs):
        """Read data from file.
        """
        method_read = getattr(phylopandas, 'read_{}'.format(schema))
        df = method_read(path, **kwargs)
        method_add = getattr(self, '_add_{}'.format(dtype))
        method_add(df)
        return self
        
    def draw_tree(self, **kwargs):
        """
        """
        if hasattr(self, 'tree') is False:
            raise Exception("TreeProject isn't away of any tree. Have you added it?")

        newick_str = self.tree.as_string(schema='newick')
        tree = toytree.tree(newick_str)
        return tree.draw(**kwargs)
        
