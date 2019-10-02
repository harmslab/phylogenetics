.. phylogenetics documentation master file, created by
   sphinx-quickstart on Wed Apr 13 13:44:41 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Phylogenetics
=============

Phylogenetics is a minimal Python API for doing phylogenetics. It manages the annoying
aspects of phylogenetics (i.e. file conversion) for you and lets you focus on exploring
and interpreting the data.

Development goals
-----------------

+ docs:
    + installation
    + specific cases
        + long branches
        + what does branch support mean?
        + sequence quality control
        + alignment quality
        + ancestors, posterior probabilities
        + how do I come up with an answerable evolutionary question?
    + design philosophy
        + why do we need another package?
        + as little low-level crap as possible (use Biopython, dendropy, etc.),
          let users interact simply via familiar csv and pandas df
+ to implement:
    + phylopandas (critical):
        + uid/csv awareness
        + tree/align integration (Zach)
    + species correction
    + datatype awareness (dna, rna, protein, codon)
    + phylobot mk 2?



Basic Example
-------------

TO DO

.. image:: _images/jlab.png
  :align: center

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
   :maxdepth: 1

   _pages/install
   _pages/project
