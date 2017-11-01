.. phylogenetics documentation master file, created by
   sphinx-quickstart on Wed Apr 13 13:44:41 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Phylogenetics
=============

Phylogenetics is a minimal Python API for doing phylogenetics. It manages the annoying
aspects of phylogenetics (i.e. file conversion) for you and lets you focus on exploring
and interpreting the data.  

Basic Example
-------------

The main object in phylogenetics is the ``TreeProject``. Start by initializing the ``TreeProject``
object, pointing to a directory where you'd like to store all the phylogenetic data/output. 

.. code-block:: python
  # Imports
  from phylogenetics import TreeProject

  # Initialize a project class
  project = TreeProject(project_dir='project')

Then, add an alignment (using any file schema you'd like). These sequences will be the
tips of your tree. You can now begin building trees and reconstructing ancestral sequences.

.. code-block:: python

  project.read_data(dtype='tips', path='alignment.fasta', schema='fasta')

  # Run PhyML to construct a phylogenetic 
  # tree by maximum likelihood.
  project.run_tree()

  # Reconstruct ancestral sequences using default settings.
  project.run_reconstruction()


``phylogenetics`` imports the ``toytree`` library to quickly plot trees in Jupyter Notebooks.
You can quickly view your tree anytime using the `draw_tree` method. 

.. code-block:: python

  project.draw_tree(width=700,
      tip_labels='id',
      tip_labels_align=True,
      use_edge_lengths=True,
      node_labels='id') 

<img src="_images/jlab.png" align="middle">

Installation
------------

Install from PyPi:

.. code-block:: bash

  pip install phylogenetics


To install a development version:

.. code-block:: bash

  git clone https://github.com/Zsailer/phylogenetics
  cd phylogenetics
  pip install -e .

Dependencies
------------

`phylogenetics` manages phylogenetics data. Currently, it doesn't do any of the phylogenetic calculations itself. For this, use external tools like:

1. PhyML_: building maximum likelihood trees.
2. PAML_: reconstructing ancestors.

.. _PhyML: http://www.atgc-montpellier.fr/phyml/
.. _PAML: http://abacus.gene.ucl.ac.uk/software/paml.html

`phylogenetics` is built on top of following python stack:

1. DendroPy_: A Python library for phylogenetic scripting, simulation, data processing and manipulation.
2. ToyTree_: A minimalist tree plotting library using toyplot graphs
3. PhyloPandas_: Pandas for phylogenetics
4. PyASR_: Ancestral Sequence Reconstruction in Python

.. _DendroPy: https://github.com/jeetsukumaran/DendroPy
.. _ToyTree: https://github.com/eaton-lab/toytree
.. _PhyloPandas: https://github.com/Zsailer/phylopandas
.. _PyASR: https://github.com/Zsailer/pyasr

Table of Contents
=================

.. toctree::
   :maxdepth: 2

   _pages/install

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
