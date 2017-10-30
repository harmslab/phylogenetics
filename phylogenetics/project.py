import os
import shutil
import subprocess
import pkg_resources
import pandas
import toytree
import phylopandas
import pyasr
import dendropy
from time import localtime, strftime

from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.PAML import codeml


def track_in_history(method):
    """Track this call in the history DataFrame"""
    def wrapper(self, *args, **kwargs):
        """"""
        # Now run method
        output = method(self, *args, **kwargs)
        
        # Store method in history.
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
        
        # Write history to file.
        history_file = os.path.join(self.project_dir, 'history.csv')
        self.history.to_csv(history_file)
        return output
    return wrapper

class TreeProject(object):
    """Main Phylogenetics Project class."""
    def __init__(self, project_dir):
        # Get current time for history.
        time = strftime("%Y-%m-%d %H:%M:%S", localtime())
        
        # Set up a project directory
        self.project_dir = project_dir
        if os.path.exists(self.project_dir):
            reply = input('Project directory already exists. Overwrite? [Y/n]')
            if reply != 'Y':
                raise Exception('Please give a new directory name.')
            else:
                shutil.rmtree(self.project_dir)
        
        # Create directory.
        os.makedirs(self.project_dir)
        
        # DataFrame for storing data.
        self.data = {'tips':None, 'ancs':None, 'tree':None}
        
        # Create a history database.
        self.history = pandas.DataFrame({'time':[time], 
            'method':['__init__'],
            'args':[project_dir], 
            'kwargs':[None]}, dtype=str)
            
        # Write history to file.
        history_file = os.path.join(self.project_dir, 'history.csv')
        self.history.to_csv(history_file)
          
    def __str__(self):
        # Get history
        history = self.history.iloc[-1]
        tips_bool = self.data['tips'] is not None
        ancs_bool = self.data['ancs'] is not None
        tree_bool = self.data['tree'] is not None
        
        # Start building history.
        info = ["TreeProject(project_dir={})\n".format(self.project_dir),
                "    last modified\t{}\n".format(history['time']),
                "    last edit\t\t{}\n".format(history['method'])]
                
        # If leafs, add stats
        info += ["    tips\t\t{}\n".format(tips_bool)]
        if tips_bool:
            leafs = self.data['tips']
            info += ["        - num of tips\t{}\n".format(len(leafs))]
    
        # If ancestors, add stats.
        info += ["    ancs\t\t{}\n".format(ancs_bool)]
        if ancs_bool:            
            ancs = self.data['ancs']
            info += ["        - num of ancs\t{}\n".format(len(ancs))]            
        
        # If tree, add stats.
        info += ["    tree\t\t{}\n".format(tree_bool)]
        return "".join(info).strip()
                
    def __repr__(self):
        return self.__str__()
    
    def save(self):
        # Save history
        history_path = os.path.join(self.project_dir, 'history.csv')
        self.history.to_csv(history_path)
        
        # Save tips data.
        if self.data['tips'] is not None:
            tips_path = os.path.join(self.project_dir, 'tips.csv')
            self.data['tips'].to_csv(tips_path)
        
        # Save ancs data.
        if self.data['ancs'] is not None:
            ancs_path = os.path.join(self.project_dir, 'ancs.csv')
            self.data['ancs'].to_csv(ancs_path)

        # Save tree data.
        if self.data['tree'] is not None:
            tree_path = os.path.join(self.project_dir, 'tree.newick')
            self.data['tree'].write(path=tree_path, schema='newick')
        return self
    
    @classmethod
    def load(cls, project_dir):
        pass
    
    @track_in_history
    def test(self, blah, blah2=None):
        return self

    def _add_tips(self, data):
        """"""
        if isinstance(data, pandas.DataFrame) == False and isinstance(data, phylopandas.DataFrame) == False:
            raise Exception('Bad datatype.')
        
        # Add unique ids
        if 'unique_id' not in data:
            unique_ids = ["tip{:07d}".format(i) for i in range(len(data))]
            col = phylopandas.Series(unique_ids, index=data.index)
            data['unique_id'] = col
        
        self.data['tips'] = data
        self.tips = data

    def _add_ancs(self, data):
        """"""
        if isinstance(data, pandas.DataFrame) == False and isinstance(data, phylopandas.DataFrame) == False:
            raise Exception('Bad datatype.')
        
        # Add unique ids
        if 'unique_id' not in data:
            unique_ids = ["anc{:07d}".format(i) for i in range(len(data))]
            col = phylopandas.Series(unique_ids, index=data.index)
            data['unique_id'] = col
        
        self.data['ancs'] = data
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
        if dtype not in ['tips', 'ancs', 'tree']:
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
    
    def write_data(self, data, path, schema='fasta', **kwargs):
        """
        """
        df = self.data[data]
        method_write = getattr(df, 'to_{}'.format(schema))
        method_write(filename=path)
        return self
    
    @track_in_history
    def run_tree(self,
        datatype='aa',
        bootstrap='0',
        model='LG',
        frequencies='e',
        ):
        """Use PhyML to build a phylogenetic tree."""
        df = self.data['tips']
        
        # Write file to disk
        fasta_file = os.path.join(self.project_dir, 'alignment.phy')
        df.to_phylip(fasta_file, sequence_col='alignment', id_col='unique_id')
        
        # Prepare options
        options = {
            'input':fasta_file,
            'datatype':datatype,
            'bootstrap':bootstrap,
            'model':model,
            'frequencies':frequencies
        }
        
        # Build command line arguments for PhyML.
        cml = PhymlCommandline(**options)
        cml_args = str(cml).split()
        output = subprocess.run(cml_args)
        
        # # Get path
        tree_file = os.path.join(self.project_dir, 'alignment.phy_phyml_tree.txt')
        tree = dendropy.Tree.get(path=tree_file, schema='newick')
        self._add_tree(tree)
        return self

    @track_in_history
    def run_reconstruction(self,
        verbose=None,
        CodonFreq=None,
        cleandata=None,
        fix_blength=None,
        NSsites=None,
        fix_omega=None,
        clock=None,
        ncatG=None,
        runmode=None,
        fix_kappa=None,
        fix_alpha=None,
        Small_Diff=None,
        method=None,
        Malpha=None,
        aaDist=None,
        RateAncestor=None,
        aaRatefile='lg',
        icode=None,
        alpha=None,
        seqtype=None,
        omega=None,
        getSE=None,
        noisy=None,
        Mgene=None,
        kappa=None,
        model=3,
        ndata=None):
        """Use PAML to build a phylogenetic tree."""
        df = self.data['tips']
        tree = self.data['tree']

        # Write file to disk
        ### NEED TO FIX SEQUENCE COL!!

        n, m = len(df), len(df['sequence'][0])
        alignment_file = 'alignment.phy'
        alignment_path = os.path.join(self.project_dir, alignment_file)
        alignment_str = df.to_fasta(sequence_col='sequence', id_col='unique_id', id_only=True)
        alignment_str = "{} {}\n".format(n,m) + alignment_str
        
        with open(alignment_path, 'w') as f:
            f.write(alignment_str)
        
        # Write tree to file
        tree_file = 'tree-to-reconstruct.newick'
        tree_path = os.path.join(self.project_dir, tree_file)
        tree.write(path=tree_path, schema='newick', 
            suppress_internal_taxon_labels=True,
            suppress_internal_node_labels=True)
        
        output_file = 'results.txt'
        output_path = os.path.join(self.project_dir, output_file)
        
        # copy model from package to project directory. 
        path_to_model = pkg_resources.resource_filename('phylogenetics', os.path.join('dat', '{}.dat'.format(aaRatefile)))
        model_file = '{}.dat'.format(aaRatefile)
        model_path = os.path.join(self.project_dir, model_file)
        shutil.copyfile(path_to_model, model_path)
    
        # Build control file.
        cml = codeml.Codeml(alignment=alignment_path, 
            tree=tree_path,
            out_file=output_path,
            working_dir=self.project_dir)
        cml.set_options(
            verbose = verbose,
            CodonFreq = CodonFreq,
            cleandata = cleandata,
            fix_blength = fix_blength,
            NSsites = NSsites,
            fix_omega = fix_omega,
            clock = clock,
            ncatG = ncatG,
            runmode = runmode,
            fix_kappa = fix_kappa,
            fix_alpha = fix_alpha,
            Small_Diff = Small_Diff,
            method = method,
            Malpha = Malpha,
            aaDist = aaDist,
            RateAncestor = RateAncestor,
            aaRatefile = model_file,
            icode = icode,
            alpha = alpha,
            seqtype = seqtype,
            omega = omega,
            getSE = getSE,
            noisy = noisy,
            Mgene = Mgene,
            kappa = kappa,
            model = model,
            ndata = ndata)
            
        # Write out control file.
        cml.ctl_file = os.path.join(self.project_dir, 'codeml_options.ctl')
        cml.write_ctl_file()     

        # Run PAML: 1. change directory, 2. run codeml, 3. parse results and 4. switch back to current dir.
        current_dir = os.getcwd()
        project_dir = os.path.join(current_dir, self.project_dir)
        os.chdir(project_dir)
        output = subprocess.run(['codeml', 'codeml_options.ctl'])
        os.chdir(current_dir)
        
        # Parse output.
        rst_file = os.path.join(self.project_dir, 'rst')
        tree_ancs, df_ancs = pyasr.read_codeml_output(rst_file)
        self._add_ancs(df_ancs)
        self._add_tree(tree_ancs)
        return self
    
    def draw_tree(self, width=200, height=500,
        tip_labels=None,
        tip_labels_color=None,
        tip_labels_align=False,
        use_edge_lengths=False,
        edge_style=None,
        edge_align_style=None,
        node_labels=None,
        node_size=None,
        node_labels_style=None,
        node_color=None,
        **kwargs):
        """Draw tree (using the toytree package.)
        """
        if hasattr(self, 'tree') is False:
            raise Exception("TreeProject isn't away of any tree. Have you added it?")

        newick_str = self.tree.as_string(schema='newick')
        tree = toytree.tree(newick_str)
        
        # Tip settings
        tree_tips_labels = tree.get_tip_labels()
        tree_ancs_labels = tree.get_node_values('name', show_root=False, show_tips=False)
        
        # Get tips
        tips = self.data['tips']
        ancs = self.data['ancs']
        
        # set tip labels to column in tip dataframe
        if type(tip_labels) is str:
            mapping = dict(zip(tree_tips_labels, tips[tip_labels]))
            tip_labels = [mapping[label] for label in tree_tips_labels]

        # set tip label colors to column in tips dataframe
        if type(tip_labels_color) is str:
            mapping = dict(zip(tree_tips_labels, tips[tip_label_colors]))
            tip_labels_color = [mapping[label] for label in tree_tips_labels]

        # # set anc labels to column in anc dataframe
        # if type(node_labels) is str:
        #     mapping = dict(zip(ancs['id'], ancs[node_labels]))
        #     node_labels = [mapping[tree_ancs_labels[i]] for i in range(len(tree_ancs_labels))]            
        # 
        # # set anc labels to column in anc dataframe
        # if type(node_size) is str:
        #     mapping = dict(zip(ancs['id'], ancs[node_size]))
        #     node_size = [mapping[tree_ancs_labels[i]] for i in range(len(tree_ancs_labels))]    
        # 
        # # set anc labels to column in anc dataframe
        # if type(node_color) is str:
        #     mapping = dict(zip(ancs['id'], ancs[node_color]))
        #     node_color = [mapping[tree_ancs_labels[i]] for i in range(len(tree_ancs_labels))]  

        options = dict(
            width=width,
            height=height,
            tip_labels=tip_labels,
            tip_labels_color=tip_labels_color,
            tip_labels_align=tip_labels_align,
            use_edge_lengths=use_edge_lengths,
            edge_style=edge_style,
            edge_align_style=edge_align_style,
            node_labels=tree_ancs_labels, #### FIX LATER
            node_size=node_size,
            node_labels_style=node_labels_style,
            node_color=node_color)

        return tree.draw(**options)
        
