# API install
from setuptools import setup

setup(name = 'phylogenetics',
    version = '0.3',
    description = 'Python API that provides simple tools for doing phylogenetics.',
    author = 'Zach Sailer',
    author_email = 'zachsailer@gmail.com',
    url = 'https://github.com/Zsailer/phylogenetics',
    download_url = 'https://github.com/Zsailer/phylogenetics/tarball/v0.3',
    packages = ['phylogenetics'],
    install_requires=[
        'dendropy',
        'biopython'
    ],
    #scripts = 'scripts/blast-download',
        #'scripts/blast-process',
        #'scripts/blast-reverse',
        #'scripts/blast-seeds',
        #'scripts/phylo-add-align',
        #'scripts/phylo-align',
        #'scripts/phylo-cdhit',
	    #'scripts/phylo-concat',
        #'scripts/phylo-edit-names',
        #'scripts/phylo-rm',
        #'scripts/phylo-tree'],
    zip_safe = False
)
