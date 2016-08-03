# Python API for doing phylogenetics

This is the master repository for the `phylogenetics` Python package. This package offers an intuitive API for analyzing, creating, and reading/writing phylogenetic data via Python. The goal of this API is to make molecular phylogenetics broadly accessible without requiring all the technical/computational expertise currently needed.

See the documentation [here](http://phylogenetics.readthedocs.io/en/latest/)

From a programming point-of-view, this API defines/standardizes a set of objects and data-structures for managing phylogenetic data. Perhaps the most frustrating aspect of doing phylogenetics (and most of bioinformatics) is converting your sequence data to various file formats used by the separate phylogenetic tools (ie. fasta, phylip, newick, nexus, etc.). One huge strength of this API is that it manages data I/O for you. It reads and writes formats seamlessly -- most of which you won't even notice as a user.

This package provides tools for:

1. Constructing and managing multiple sequence alignments (using MSAProbs).
2. Constructing and analyzing maximum likelihood phylogenetic trees (using PhyML and DendroPy).
3. Reconstructing ancestral sequences via ASR (via PAML).
4. Reading and Writing to various file formats to standing Python objects.

## The Basics

The main entry-point to the API is the `Project` object, which acts as a container object for all pieces of a
phylogenetics project. Below is an example of how to use the `Project` class.

```python
from phylogenetics.project import Project

# Initialize a project class
project = Project()
```
Now we need to start adding data to this empty phylogenetics object. We can do this by downloading sequences
from BLAST database.
```python
# List some accession IDS to download from BLAST database.
email = ""
accessions = [
    "AGH62057"
    "NP_004553",
    "AHW56551",
    "BAA25751",
    "ABN46990"
]

# Download from database!
project.download(accession, email)
```
Alternatively, you may have sequence data in a file somewhere. In which case, you can load them into your project by reading them:
```python
project.Read.fasta(fname="sequence.fasta")
```
The `Read` module read many different formats.

At this point, the sequences and metadata for the accession IDs listed are now
contained by `project` in an object called `HomologSet`.

We can now align the `HomologSet`, build a tree, and reconstruct the it's ancestors by following methods:
```python
project.align()
project.tree()
project.reconstruct()
```
Doing a full phylogenetics project is as simple as that!

To see more, check out some Jupyter Notebook [examples](https://github.com/Zsailer/phylogenetics/tree/master/examples).

## Coming Soon

A few features on the immediate roadmap for this project are:

1. Plotting modules for analyzing pieces of phylogenetic projects on the fly (via `matplotlib`)
2. Interactive widgets for manipulating data on the fly (via `ipywidgets`)
3. Light-weight treeviewer to quickly visualize trees.

## Installation

### User
This package can be installed via PyPI:
```
pip install phylogenetics
```
Note, however, that development might be changing rapidly, and the builds on PyPI might
fall behind quickly. To keep up with development, clone this repo locally:
```
git clone https://github.com/Zsailer/phylogenetics
```
Navigate to this directory, and install this python package with
```
python setup.py install
```

### Development
If you'd like to work on developing `phylogenetics`, install with:
```
pip install -e .
```

**NOTE:** Many of modules in this API require other software packages to work.

## Dependencies

`phylogenetics` is an API for managing phylogenetics data. It does not, on its own, run any of the
calculations. Instead, it offers wrappers around current phylogenetic tools that must be
downloaded and built separately. The [wiki](https://github.com/Zsailer/phylogenetics/wiki/Setting-up-Phylogenetics-package) provides some documentation on how to install many of
these dependencies.

External dependencies include:

1. `cdhit` - clustering redundant sequences
2. `msaprobs` - multiple sequence alignment software
3. `phyml` - building maximum likelihood trees
4. `paml` - reconstructing ancestors

Python dependencies includes:

1. `biopython` (for Entrez calls).
2. `dendropy`  -- and AWESOME API for working with Tree data-structures.
3. `requests` -- making requests to web service API's.

## Credits

Most of this work was originally inspired [Mike Harms'](https://github.com/harmsm) `phylo_tools` repo. The `Reconstruction` modules of this package were inspired by [Victor Hanson-Smith's](https://github.com/vhsvhs) Python package `lazarus`. While Lazarus is not used by this API, much the work for doing ASR in Python was pioneered by this process.
