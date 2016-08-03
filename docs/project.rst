Project Class
=============

The ``project`` class provides a single, persisent object that walks you through
a phylogenetics pipeline. As you construct your project, new methods and attributes emerge
to analyze your sequence data.

This object also acts as a container for all your data, and saves its current state
at each step. (To see what a saved project looks like, checkout the the `Metadata page`_.)
The ```phylogenetics`` API, then, makes it easy to walk through your data/analysis.

.. _Metadata page:

Contents of this page:

* `Project design`_
* `Basic Usage`_


Project design
--------------

The diagram below shows the architecture of the ``Project`` class and its API.

.. image:: _images/project-fig.png

Each object inside the Project class attaches a set of useful methods for analyzing
each type of data. Each piece of a phylogenetics package lives shallow in the project
class, and can be accessed through this single object's API. The modular design
of this package also enables you to build each piece of a project separately. You
can always add it to a new/old project later!

Basic Usage
-----------
There a few ways in which you can construct a phylogenetics Project.

- Start fresh
    1. initialize a project object::

        project = phylogenetics.Project()

    2. download a list of accession ids::

        project.download(accessions_list)

    3. run pipeline::

        project.cluster()
        project.align()
        project.tree()
        project.reconstruct()

- Start with existing `phylogenetics` objects/data

    1. Initialize project class with other objects given as arguments::

        project = phylogenetics.Project(
            HomologSet,
            Alignment,
        )

    2. Continue through pipeline

    3. Or add separate objects later.::

        project = phylogenetics.add(Tree)

- Load existing phylogenetics file from file.

    Each action in the phylogenetics package saves a copy to disk by default.::

        project.load("project-#####.pickle")

- Read phylogenetic data from files into project class.

    1. Initialize a project object::

        project = phylogenetics.Project()

    2. Read files containing phylogenetic information using the `files` method.::

        project.files(
            HomologSet="sequences.fasta",
            Alignment="alignment.fasta",
            Tree="tree.newick",
        )

    3. Or read in files individually.::

        project.Read.fasta("sequences.fasta")
