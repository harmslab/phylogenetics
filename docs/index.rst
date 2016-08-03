.. phylogenetics documentation master file, created by
   sphinx-quickstart on Wed Apr 13 13:44:41 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Phylogenetics package
=====================

This package was designed for two purposes:

1. Create a single, persistent python object for managing phylogenetic data
2. Simple, one-stop python API for running through a complete phylogenetic pipeline.

Because of the API nature of this package, your phylogenetics analysis can be
easily run in a `Jupyter Notebook`_, making your analysis clean, clear, and reproducible.

Try it out!
-----------

* View static `examples`_ or
* `demo`_ the package using `Binder`_ notebooks.

.. _Jupyter Notebook: http://jupyter.org/
.. _examples: http://nbviewer.jupyter.org/github/Zsailer/phylogenetics/tree/master/examples/
.. _demo: http://mybinder.org/repo/Zsailer/phylogenetics
.. _Binder: http://mybinder.org/

What does it do?
----------------

Using external tools (see how to install them `here`_), this package makes the
following tasks easy to do programmatically.

* BLAST a set of sequences
* Download taxonomic data
* Align a set of sequences
* Construct a phylogenetic tree from sequences
* Reconstruct ancestors in a phylogenetics tree
* Read and write many file formats in phylgenetics
* Manipulate, combine, subtree, etc. phylogenetic trees.

.. _here: https://github.com/Zsailer/phylogenetics/wiki/Setting-up-Phylogenetics-package

Table of Contents
=================

.. toctree::
   :maxdepth: 2

   install
   project
   reading
   writing
   metadata

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
