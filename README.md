# Python API with tools for doing phylogenetics

This is the master repository for the `phylogenetics` Python package. This package includes a lightweight, simple-to-use API for managing and processing phylogenetic data. Many of the modules included in this package originated from Dr. [Mike Harms'](https://github.com/harmsm) `phylo_tools` and have been converted to API's. Read the Wiki to learn more about the internal structure. 

The foundation of this API are the `Homolog` and `HomologSet` objects. These objects offer a simple datastructure that manages the metadata for a set of sequences in a phylogenetics/reconstruction project. These objects are easily queried, updated, and saved into many formats (i.e. fasta, csv, phylip, pickle, and json). 

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

Clone this repo locally:

```
git clone https://github.com/Zsailer/phylogenetics
```

Navigate to this directory, and install this python package with

```
python setup.py install
```

**NOTE:** Many of modules in this API require

## Setting up for Development

Git must be installed to clone and contribute to this project


1. Fork this repository on Github
2. Clone that repository locally
```
git clone https://github.com/Zsailer/phylogenetics
```
3. Navigate to this directory, and install (softly) this python package with
```
cd phylogenetics
python setup.py develop
```
4. Add another remote link to the master version, call it `upstream`.
```
git remote add upstream
```
5. Start a branch locally from local master
```
git checkout -B <branch-name>
```
6. Make changes and commit to that branch.
```
git commit -a -m "<commit message>"
```
7. Push to your fork on github (which you called `upstream`).
```
git push upstream <branch-name>
```
8. Pull request the branch on Github into this master repository on Github.

