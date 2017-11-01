# Python API for managing a phylogenetics project

[![Documentation Status](http://readthedocs.org/projects/phylogenetics/badge/?version=latest)](http://phylogenetics.readthedocs.io/en/latest/?badge=latest)


Phylogenetics is a minimal Python API for doing phylogenetics. It manages the annoying aspects of phylogenetics (i.e. file conversion) for you and lets you focus on exploring and interpreting the data.  

**Note** The `phylogenetics` package has been completely rewritten between v0.3 and v0.4. v0.3 is now deprecated and no longer maintained. v0.4 has significantly simplified the API. I hope you enjoy it.

## Basic Example

The main object in phylogenetics is the `TreeProject`. Start by initializing the `TreeProject`
object, pointing to a directory where you'd like to store all the phylogenetic data/output. 

```python
# Imports
from phylogenetics import TreeProject

# Initialize a project class
project = TreeProject(project_dir='project')
```

Then, add an alignment (using any file schema you'd like). These sequences will be the
tips of your tree. You can now begin building trees and reconstructing ancestral sequences.

```python
project.read_data(dtype='tips', path='alignment.fasta', schema='fasta')

# Run PhyML to construct a phylogenetic 
# tree by maximum likelihood.
project.run_tree()

# Reconstruct ancestral sequences using default settings.
project.run_reconstruction()
```

`phylogenetics` imports the `toytree` library to quickly plot trees in Jupyter Notebooks.
You can quickly view your tree anytime using the `draw_tree` method. 
```python
project.draw_tree(width=700,
    tip_labels='id',
    tip_labels_align=True,
    use_edge_lengths=True,
    node_labels='id')
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

3. [PhyML](http://www.atgc-montpellier.fr/phyml/) - building maximum likelihood trees.
4. [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) - reconstructing ancestors.

`phylogenetics` is built on top of following python stack:

1. Pandas 
2. Biopython
3. DendroPy
4. ToyTree
5. PhyloPandas
6. PyASR
