# API install
from setuptools import setup

setup(name = 'phylogenetics',
    version = '1.0',
    description = 'Human readable phylogenetics in Python.',
    author = 'Zach Sailer',
    author_email = 'zachsailer@gmail.com',
    url = 'https://github.com/Zsailer/phylogenetics',
    packages = ['phylogenetics'],
    install_requires=[
        'dendropy',
        'biopython'],
    zip_safe = False
)
