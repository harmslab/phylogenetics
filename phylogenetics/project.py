import os
import pickle
import subprocess
import pandas as pd
import phylopandas as ph
import phylogenetics.tools

import pyasr
import Bio

from .history import track_in_history

class PhylogeneticsProject(object):
    """A lightweight Python object that populates a PhyloPandas DataFrame.

    Parameters
    ----------
    project_dir : str
        the directory to store the phylogenetic data.

    overwrite : bool (default: False)
        allow overwriting a project that already exists in project_dir location.
    """
    @track_in_history
    def __init__(self, project_dir, overwrite=False):

        # Set up a project directory
        if os.path.exists(project_dir) and overwrite is False:
            raise Exception("Project already exists! Use `PhylogeneticsProject.load_pickle` or delete the project.")
        elif not os.path.exists(project_dir):
            os.makedirs(project_dir)

        # Define columns for project.
        columns = [
            'uid',
            'description',
            'id',
            'label',
            'sequence',
            'type',
            'parent',
            'branch_length'
        ]

        self.project_dir = project_dir
        self.data = pd.DataFrame(columns=columns)
        self.history = list()

    @staticmethod
    def load_pickle(filename):
        """"""
        with open(filename, 'rb') as f:
            self = pickle.load(f)
        return self

    def to_pickle(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


    @track_in_history
    def read_data(self, path, schema, **kwargs):
        """
        """
        method_name = 'read_{}'.format(schema)
        try:
            method = getattr(self.data.phylo, method_name)
        except:
            method = getattr(self.data, method_name)
        self.data = method(path, **kwargs)

    @track_in_history
    def compute_blast(self, **kwargs):

        print("Running blast...")

        new_df = phylogenetics.tools.blast.ncbi_blast(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + "/0_blast_output.csv")

        print("Done.")

    @track_in_history
    def compute_clusters(self, **kwargs):

        print("Computing clusters...")

        new_df = phylogenetics.tools.cdhit.run(self.data)
        self.data = new_df

        new_df.to_csv(self.project_dir + "/1_clustering_output.csv", **kwargs)

        print("Done.")


    @track_in_history
    def compute_alignment(self, **kwargs):

        print("Computing alignment...")

        new_df = phylogenetics.tools.align.run(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + "/2_alignment_output.csv")

        print("Done.")

    @track_in_history
    def compute_gblocks(self, **kwargs):

        print("Computing gblocks...")

        new_df = phylogenetics.tools.gblocks.run(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + "/3_gblocks_alignment_output.csv")

        print("Done.")


    @track_in_history
    def compute_tree(self, **kwargs):

        print("Computing tree...")

        new_df = phylogenetics.tools.tree.run(self.data, self.project_dir)
        self.data = new_df

        new_df.to_csv(self.project_dir + "/4_tree_output.csv", **kwargs)

        print("Done.")

    @track_in_history
    def compute_reconstruction(self, **kwargs):

        print("Computing reconstruction...")

        new_df = phylogenetics.tools.ancestral_reconstruction.run(self.data, self.project_dir)
        self.data = new_df

        new_df.to_csv(self.project_dir + "/5_reconstruction_output.csv", **kwargs)

        print("Done.")

    @track_in_history
    def fetch_previous_project(self, old_data, ancestors=None, **kwargs):
        """
        Read csv(s) from previous project and using pandas.read_csv function
        old_data: csv from previous project
        """

        self.data = pd.read_csv(old_data)
        new_df.to_csv(self.project_dir + "/old_data.csv")
