"""Phylogenetics package entry point.

`Project` class is a versatile object that manages all data for a phylogenetics
project.

Example:

>>> Project.Alignment
            .Tree
            .HomologSet
            .AncestorSet


"""
import os, pickle, datetime

# import objects to bind to Project class
from phylogenetics.homologs import Homolog, HomologSet
from phylogenetics.alignment import Alignment
from phylogenetics.tree import Tree
from phylogenetics.ancestors import Ancestor, AncestorSet
from phylogenetics.reconstruction import Reconstruction
from phylogenetics.dataio import projectio

# imports for running external tools.
from .exttools import (cdhit,
                        msaprobs,
                        phyml,
                        paml)


class Project(object):
    """Container object for managing all data for a phylogenetics project.

    Optional arguments can be passed into the Project class. These must be
    objects from the `phylogenetics` package (i.e. HomologSet, Alignment, Tree,
    Reconstruction, AncestorSet, etc.)


    There a few ways in which you can construct a project object.

    1. Start fresh
        - initialize a project object
            >>> project = phylogenetics.Project()
        - download a list of accession ids
            >>> project.download(accessions_list)
        - run pipeline,
            >>> project.cluster()
            >>> project.align()
            >>> project.tree()
            >>> project.reconstruct()

    2. Start with existing `phylogenetics` objects/data
        - Initialize project class with other objects given as arguments
            >>> project = phylogenetics.Project(
                            HomologSet,
                            Alignment,
                        )
        - Continue through pipeline
        - Or add separate objects later.
            >>> project = phylogenetics.add(Tree)

    3. Load existing phylogenetics file from file.
        - Each action in the phylogenetics package saves a copy to disk by default.
            >>> project.load("project-#####.pickle")

    4. Read phylogenetic data from files into project class.
        - Initialize a project object
            >>> project = phylogenetics.Project()
        - Read files containing phylogenetic information using the `files` method.
            >>> project.files(
                HomologSet="sequences.fasta",
                Alignment="alignment.fasta",
                Tree="tree.newick",
            )
        - Or read in files individually.
            >>> project.Read.fasta("sequences.fasta")

    """
    def __init__(self, *args):
        # components
        self._components = {}
        # Bind Reading module to class
        self.Read = projectio.Read(self)
        self.Write = projectio.Write(self)
        # Add any objects that were given to Project
        for a in args:
            self.add(a)

    @classmethod
    def load(cls, path):
        """ Load a project from pickle file. """
        project = cls()
        project.Read.pickle(fname=path)
        return project

    def save(self, fname="project-%s.pickle" % \
        datetime.date.today().isoformat()):
        """ Save Project to path. """
        self.Write.pickle(fname=fname)

    def files(self, **files):
        """Load files into project object.

        **files are keyword arguments.
        """
        # Object types accessible to Project object
        object_types = {
            "HomologSet" : HomologSet,
            "Alignment" : Alignment,
            "Tree" : Tree,
            "AncestorSet" : AncestorSet,
            "Reconstruction" : Reconstruction,
        }

        # File extensions used by project object
        ext_types = {
            "fasta" : "fasta",
            "json" : "json",
            "nwk" : "newick",
            "nxs" : "nexus",
            "pickle" : "pickle",
            "rst" : "rst",
            "csv" : "csv",
            "xml" : "entrez_xml",
            "phy" : "phylip"
        }

        for key, fname in files.items():
            # Find file-type
            filename, extension = os.path.splitext(fname)
            filetype = ext_types[extension[1:]]
            # Initialize the object that this file populates
            new_object = object_types[key]()
            # Get the reading method from new object
            read_method = getattr(new_object.Read, filetype)
            # Read into that object
            read_method(fname=fname)
            # Add object to project object.
            self.add(new_object)

    def add(self, item):
        """Add data to project.
        """
        # possible objects to add
        items = {
            HomologSet: self._add_HomologSet,
            Alignment: self._add_Alignment,
            Tree: self._add_Tree,
            Reconstruction: self._add_Reconstruction,
            AncestorSet: self._add_AncestorSet
        }

        # Find item type in set of possible items
        adding_method = items[item.__class__]

        # Add that item to project
        adding_method(item)

    def download(self, ids, email):
        """ Download a set of Homologs"""
        hs = HomologSet()
        self._add_HomologSet( hs )
        self.HomologSet.download(ids, email)

    def _add_HomologSet(self, HomologSet):
        """Add HomologSet Set to PhylogeneticsProject object."""
        self._components["HomologSet"] = HomologSet
        # Set the HomologSet object
        self.HomologSet = HomologSet
        # Expose the align method of this object to user
        setattr(self, "cluster", self._cluster)
        setattr(self, "align", self._align)

    def _add_Alignment(self, Alignment):
        """Add Alignment to PhylogeneticsProject object."""
        self._components["Alignment"] = Alignment
        # Set the Alignment object
        self.Alignment = Alignment
        # Expose the tree methods of this project object
        setattr(self, "tree", self._tree)

    def _add_Tree(self, Tree):
        """Add Tree to PhylogeneticsProject object."""
        self._components["Alignment"] = Tree
        # Set the Tree object of project
        self.Tree = Tree
        # Expose the reconstruction methods of this project object
        setattr(self, "reconstruct", self._reconstruct)

    def _add_Reconstruction(self, Reconstruction):
        """Add Reconstruction to PhylogeneticsProject object."""
        self._components["Reconstruction"] = Reconstruction
        self.Reconstruction = Reconstruction

    def _add_AncestorSet(self, AncestorSet):
        """Add a AncestorSet object to PhylogeneticsProject object."""
        self._components["AncestorSet"] = AncestorSet
        self.AncestorSet = AncestorSet

    def _cluster(
        self,
        redund_cutoff=0.99,
        tmp_file_suffix="oB_cdhit",
        word_size=5,
        cores=1,
        keep_tmp=False,
        accession=(),
        positive=(),
        negative=("putative","hypothetical", "unnamed", "possible", "predicted",
                    "unknown", "uncharacterized","mutant", "isoform"),
        ):
        """Remove redundant sequences from HomologSet based on some sequence redundancy
        cutoff threshold.

        This method moves the original HomologSet to a new object inside project
        called `HomologSet_original`. Clustering always happens on this set. The
        results are set as the new HomologSet object.
        """
        # If the original set is old, use it for clustering
        if hasattr(self, "HomologSet_original"):
            HomologSet_to_cluster = self.HomologSet_original
        # Else HomologSet is the original set.
        else:
            HomologSet_to_cluster = self.HomologSet
            setattr(self, "HomologSet_original", self.HomologSet)

        # Cluster original set
        new_HomologSet = HomologSet_to_cluster.cluster(
            redund_cutoff=redund_cutoff,
            tmp_file_suffix=tmp_file_suffix,
            word_size=word_size,
            cores=cores,
            keep_tmp=keep_tmp,
            accession=accession,
            positive=positive,
            negative=negative,
            inplace=False
        )

        # Set the resulting subset HomologSet as the new HomologSet.
        self.HomologSet = new_HomologSet
        self.save()

    def _align(self,
        fname="alignment.fasta",
        rm_tmp=True,
        quiet=False,
        cores=2
        ):
        """ Multiple sequence alignment of the HomologSet.

            Currently, only option is to use MSAProbs.
        """
        # Write out alignment file
        self.HomologSet.Write.fasta(fname="alignment.fasta")

        # Run the alignment with MSAProbs
        output_fname = msaprobs.run(fasta_fname="alignment", cores=cores, rm_tmp=rm_tmp)

        # Attach an alignment object to HomologSet
        self._add_Alignment(Alignment(self.HomologSet))

        # Read alignment from output fasta and manage with Alignment object
        self.Alignment.Read.fasta(fname=output_fname)

        # Let us know when finished
        if quiet is False:
            print("Alignment finished.")

        # Remove fasta file.
        if rm_tmp:
            os.remove(output_fname)
        self.save()

    def _tree(self, **kwargs):
        """ Compute the maximum likelihood phylogenetic tree from
            aligned dataset.

        """
        # Write the HomologSet out as a phylip.
        self.Alignment.Write.phylip(fname="ml-tree.phy")

        # Run phyml and parse results.
        tree, stats = phyml.run("ml-tree", **kwargs)

        # Add Tree object to HomologSet
        self._add_Tree(Tree(self.HomologSet, tree, stats=stats))
        self.save()


    def _reconstruct(self):
        """ Resurrect Ancestors on Tree.
        """
        # Bind Ancestor Objects to each internal node.
        ancestors = []
        for node in self.Tree.Dendropy.internal_nodes():
            id = node.label
            ancestors.append( Ancestor(id, self.Tree))

        # Bind an AncestorSet object to HomologSet
        self._add_AncestorSet( AncestorSet(self.Tree, ancestors=ancestors) )
        self.AncestorSet._nodes_to_ancestor()

        seqfile = "asr-alignment.fasta"
        outfile = "asr-output"
        treefile = "asr-tree.nwk"

        # Prepare input files for PAML
        self.Tree.Dendropy.write(path=treefile, schema="newick", suppress_internal_node_labels=True)
        self.Alignment.Write.fasta(fname=seqfile)

        # Construct a paml job
        paml_job = paml.CodeML(
            seqfile=seqfile,
            outfile=outfile,
            treefile=treefile,
            fix_alpha=True,
            alpha=self.Tree.stats["Gamma shape parameter"],
        )

        reconstruction = Reconstruction(self.Alignment, self.Tree, self.AncestorSet, paml_job)
        self._add_Reconstruction( reconstruction )

        # Run the PAML job
        self.Reconstruction.paml_job.run()

        # Read the paml output and bind data to tree
        self.AncestorSet.Read.rst(fname="rst")

        # Infer gaps.
        self.Reconstruction.infer_gaps()
        self.save()
