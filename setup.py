# API install
from distutils.core import setup

setup(name = 'phylogenetics',
    version = '0.1',
    description = 'Python API that provides simple tools for doing phylogenetics.',
    author = 'Zach Sailer',
    author_email = 'zachsailer@gmail.com',
    url = 'https://github.com/Zsailer/phylogenetics',
    download_url = 'https://github.com/Zsailer/phylogenetics/tarball/0.1',
    packages = ['phylogenetics'],
    scripts = ['scripts/blast-seeds',
            'scripts/blast-reverse',
            'scripts/blast-process',
            'scripts/blast-download',
            'scripts/phylo-edit-names',
            'scripts/phylo-add-align',
            'scripts/phylo-align',
            'scripts/phylo-cdhit',
            'scripts/phylo-tree'
    ],
    install_requires=[
        'numpy',
    ],
    zip_safe = False
)
