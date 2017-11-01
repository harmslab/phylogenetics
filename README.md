# Python API for managing a phylogenetics project

Phylogenetics is a minimal Python API for doing phylogenetics. It manages the annoying aspects of phylogenetics (i.e. file conversion) for you and lets you focus on exploring and interpreting the data.  

**Note** The `phylogenetics` package has been completely rewritten between v0.3 and v0.4. v0.3 is now deprecated and no longer maintained. v0.4 has significantly simplified the API. I hope you enjoy it.

## The Basics

The main entry-point to the API is the `Project` object, which acts as a container object for all pieces of a
phylogenetics project. Below is an example of how to use the `Project` class.

```python
# Imports
from phylogenetics import TreeProject

# Initialize a project class
project = TreeProject(project_dir='project')
project.read_data(dtype='tips', path='alignment.fasta', schema='fasta')

# Run PhyML to construct a phylogenetic 
# tree by maximum likelihood.
project.run_tree()

# Reconstruct ancestral sequences.
project.run_reconstruction(alpha=1.2)
```

<img src="docs/_images/jlab.png" align="middle">

## Installation

Install from PyPi:
```
pip install phylogenetics
```

To install a development version:
```
git clone https://github.com/Zsailer/phylogenetics
cd phylogenetics
pip install -e .
```

## Dependencies

`phylogenetics` manages phylogenetics data. Currently, it doesn't do any of the phylogenetic calculations itself. For this, use external tools like:

1. `cdhit` - clustering redundant sequences
2. `msaprobs` - multiple sequence alignment software
3. `phyml` - building maximum likelihood trees
4. `paml` - reconstructing ancestors

`phylogenetics` is built on top of following python stack:

1. Pandas 
2. Biopython
3. DendroPy
4. ToyTree
5. PhyloPandas
