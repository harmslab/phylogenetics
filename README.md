# Python API for managing a phylogenetics project

This package provide a minimalist 

## The Basics

The main entry-point to the API is the `Project` object, which acts as a container object for all pieces of a
phylogenetics project. Below is an example of how to use the `Project` class.

```python
# Imports
from phylogenetics import TreeProject

# Initialize a project class
project = TreeProject()
```

Add 

## Installation

`phylogenetics` has been rewritten from scratch. The version on PyPi is outdated
and will not work with the new version. For now, install from source. 
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

Python dependencies includes:

1. Pandas 
2. Biopython
3. DendroPy
4. ToyTree
5. PhyloPandas
