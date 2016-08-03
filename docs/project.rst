Project Class
=============

The ``project`` class provides a single, persisent object that walks you through
a phylogenetics pipeline. As you construct your project, new methods and attributes emerge
for navigating your data.

This object acts as a container for all your data, and saves its current state
at each step.

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
