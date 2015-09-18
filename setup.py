# API install
from setuptools import setup

setup(name = 'phylogenetics',
    version = '0.1',
    description = 'Python API that provides simple tools for doing phylogenetics.',
    author = 'Zach Sailer',
    author_email = 'zachsailer@gmail.com',
    url = 'https://github.com/Zsailer/phylogenetics',
    download_url = 'https://github.com/Zsailer/phylogenetics/tarball/0.1',
    packages = ['phylogenetics'],
    scripts = ['scripts/blast_seeds.py',
            'scripts/edit_names.py',
            'scripts/reverse_blast.py',
            'scripts/xml2homologs.py'
    ],
    install_requires=[
        'numpy',
    ],
    zip_safe = False
)
