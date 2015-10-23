# Python API and command-line for doing phylogenetics

Test out the API in notebooks -- click on the badge:

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/Zsailer/phylogenetics)

This is the master repository for the `phylogenetics` Python package. This package includes a lightweight, simple-to-use API for managing and processing phylogenetic data. Many of the modules included in this package originated from Dr. [Mike Harms'](https://github.com/harmsm) `phylo_tools` and have been converted to API's. Read the Wiki to learn more about the internal structure.

This package also contains a series of command-line tools that are installed globally and do many vanilla tasks common to phylogenetics projects. While I recommend using the API in an interactive environment, like IPython or Jupyter notebooks, it is sometimes faster to use these scripts. Note that there are example notebooks included in the `examples` folder.

The foundation of this API are the `Homolog` and `HomologSet` objects. These objects offer simple data-structures that manages the metadata for a set of sequences in a phylogenetics/reconstruction project. These objects are easily queried, updated, and saved into many formats (i.e. fasta, csv, phylip, pickle, and json).

An metadata inside a homolog (__dict__ attribute) might look might look like this:

```python
>>> print(homolog.__dict__)

{
    "id" : "XX00000000",
    "species" : "S100A5",
    "organism" : "human",
    "length" : 180,
    "sequence" : "ASKFAFELGSKADGKASEKA...",
    "latest_align" : "----ASK---F-AFELG--SKA---DGKASEK-A---...",
    "align_1" : "--------------ASK-------F-AFELG-----SKA--DFLASEK--A--...",

    .
    .
    .
}

>>> homolog.fasta()

>XX00000000
ASKFAFELGSKADGKASEKA...


>>> homolog.write("homolog.fasta", format="fasta")  # writes to file "homolog.fasta"
```

## Installation

This package is on Pypi, so this works for installation:
```
pip install phylogenetics
```
however, it will likely be trailing development as I move quickly on improving this tool.

I suggest installing from source until I get a fully stable version up and running.

Clone this repo locally:
```
git clone https://github.com/Zsailer/phylogenetics
```

Navigate to this directory, and install this python package with

```
python setup.py install
```

If you'd like to work in development mode,

```
pip install -e .
```

**NOTE:** Many of modules in this API require other software packages to work.
