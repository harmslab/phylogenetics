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
        overwrite a project that already exists in project_dir location.
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
        Read sequences into project object.

        Parameters
        ----------
        path : str
            path to file containing initial seed sequence(s).
        schema : str
            type of file to read in. Most common is "fasta".
        """
        method_name = 'read_{}'.format(schema)
        try:
            method = getattr(self.data.phylo, method_name)
        except:
            method = getattr(self.data, method_name)
        self.data = method(path, **kwargs)

    @track_in_history
    def compute_blast(self, **kwargs):

        """
        Run a BLAST search using biopython BIO module.

        See https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide
        for an exhaustive description of BLASTing parameters
        and https://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-module.html
        for Biopython's NCBIWWW qblast module description.

        Parameters
        ----------
        db : str (default: "nr")
            NCBI database to search. Default non-redundant.
            Other options:
        hitlist_size : int (default: 100)
            number of sequences to return from BLAST search.
        entrez_query : str (default: '(none)')
            limit search to specific species/groups using entrez names and/or tax ids.
        e_value_cutoff : float (default: 0.01)
            expected number of chance matches in a random model.
        gapcosts : tuple (default (11,1))
            cost to create and extend a gap in an alignment.
            default 11 existence penalty, 1 extension penalty.
        """


        print("Running blast...")
        output = "/0_blast_output.csv"

        new_df = phylogenetics.tools.blast.ncbi_blast(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + output)

        print("Writing BLAST results to " + self.project_dir + output + "...Done.\n")

    @track_in_history
    def compute_clusters(self, **kwargs):

        """
        Cluster sequences by similarity cutoff using cdhit. http://weizhongli-lab.org/cd-hit/

        Sequences are sorted from longest to shortest before running cd-hit,
        ensuring that the larger sequences always end up in the final output.

        Parameters
        ----------
        cutoff : float
            percent sequence identity cutoff (between 0 and 1.0) at which to
                place two sequences in the same cluster
        """

        print("Computing clusters...")
        output = "/1_clustering_output.csv"

        new_df = phylogenetics.tools.cdhit.run(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + output)

        print("Writing clustering results to " + self.project_dir + output + "...Done.\n")

    @track_in_history
    def compute_alignment(self, **kwargs):

        """
        Align sequences using muscle (default) or msaprobs via Biopython.

        Biopython: https://biopython.org/DIST/docs/api/Bio.Align.Applications._Muscle.MuscleCommandline-class.html
        muscle source: http://www.drive5.com/muscle/

        Parameters
        ----------
        program : str (default="muscle")
            alignment program to use. Other current option is msaprobs.
        """

        print("Computing alignment...")
        output = "/2_alignment_output.csv"

        new_df = phylogenetics.tools.align.run(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + output)

        print("Writing alignment results to " + self.project_dir + output + "...Done.\n")

    @track_in_history
    def compute_gblocks(self, **kwargs):

        """Run Gblocks on a phylopandas dataframe.

        Full docs for software are here:
        http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html#Installation

        As of this writing (August 2019), Gblocks can be installed using
        bioconda (https://bioconda.github.io/).

        sequence_type: type of sequence
                       default: protein
                       allowed: ("protein","dna","codon")
                       gblock flag: -t
        min_conserved: minimum number of sequences for conserved position
                       default: num_seqs*0.5 + 1
                       allowed: (num_seqs*0.5 + 1,num_seqs)
                       gblock flag: -b1
        min_flank: minimum number of sequences for a flank position
                       default: 0.85*num_seqs
                       allowed: (min_conserved,num_seqs)
                       gblock flag: -b2
        max_contig_noncon: maximum number of contiguous non-conserved positions
                       default: 8
                       allowed: any integer > 0
                       gblock flag: -b3
        min_block_length: minimum length of a block
                       default: 10
                       allowed: any integer > 2
                       gblock flag: -b4
        allowed_gap: allowed gap positions
                       default: "no"
                       allowed: ("no","half","all")
                       gblock flag: -b5
        use_sim_matrix: use similarity matrix
                       default: True
                       allowed: True,False (True only allowed for proteins)
                       gblock flag: -b6
        min_initial_block: minimum length of an initial block
                       default: min_block_length
                       allowed: any integer > 2
                       gblock flag: -b0
        keep_tmp: keep temporary files and write out reports
                       default: False
                       allowed: False,True
        """

        print("Computing gblocks...")
        output = "/3_gblocks_alignment_output.csv"

        new_df = phylogenetics.tools.gblocks.run(self.data, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + output)

        print("Writing gblocks results to " + self.project_dir + output + "...Done.\n")

    @track_in_history
    def compute_tree(self, **kwargs):

        """
        Compute tree using PhyML.

        http://www.atgc-montpellier.fr/phyml/
        Handy user manual for PhyML: https://github.com/stephaneguindon/phyml/blob/master/doc/phyml-manual.pdf

        Parameters
        ----------
        sequence_col : str (default='sequence')

        datatype : str (default='aa')
            type of data. 'aa' for amino acids, 'dna' for DNA.
        bootstrap: str (default='-1')
            approximate likelihood ratio test - returns aLRT statistics.
        model: str (default='LG')
            substitution model used. There are many to choose from -
            see PhyML site above for other options.
        frequencies: str (default='e')
            character (amino acid or DNA) frequencies.
            Estimated by counting occurences in alignment.
        **kwargs
        """

        print("Computing tree...")
        output = "/4_tree_output.csv"

        new_df = phylogenetics.tools.tree.run(self.data, self.project_dir, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + output)

        print("Writing tree results to " + self.project_dir + output + "...Done.\n")

    @track_in_history
    def compute_reconstruction(self, **kwargs):

        """
        Run ancestral sequence reconstruction using PAML.

        http://abacus.gene.ucl.ac.uk/software/paml.html
        User guide: http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf

        Parameters
        ----------
        altall_cutoff : float (default=0.2)
            cutoff value (between 0-1) for Alt-All alternate reconstructions.
            See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6095102/
        aaRatefile : str (default='lg')
            Amino acid substitution rate matrix.
            See user guide above for other exhaustive options.
        """


        print("Computing reconstruction...")
        output = "/5_reconstruction_output.csv"

        new_df = phylogenetics.tools.ancestral_reconstruction.run(self.data, self.project_dir, **kwargs)
        self.data = new_df

        new_df.to_csv(self.project_dir + output)

        print("Writing reconstruction results to " + self.project_dir + output + "...Done.\n")

    @track_in_history
    def fetch_previous_project(self, old_data, ancestors=None, **kwargs):
        """
        Read csv(s) from previous project and using pandas.read_csv function
        old_data: csv from previous project
        """

        self.data = pd.read_csv(old_data)
        new_df.to_csv(self.project_dir + "/old_data.csv")
