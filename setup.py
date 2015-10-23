# API install
from setuptools import setup

setup(name = 'phylogenetics',
    version = '0.2.0a1',
    description = 'Python API that provides simple tools for doing phylogenetics.',
    author = 'Zach Sailer',
    author_email = 'zachsailer@gmail.com',
    url = 'https://github.com/Zsailer/phylogenetics',
    download_url = 'https://github.com/Zsailer/phylogenetics/tarball/0.1',
    packages = ['phylogenetics'],
    scripts = ['scripts/blast-download',
        'scripts/blast-process',
        'scripts/blast-reverse',
        'scripts/blast-seeds',
        'scripts/phylo-add-align',
        'scripts/phylo-align',
        'scripts/phylo-cdhit',
        'scripts/phylo-edit-names',
        'scripts/phylo-tree'],
    zip_safe = False
)
