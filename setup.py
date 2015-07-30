# API install
from setuptools import setup

setup(name='phylogenetics',
    version='0.1',
    description='Python API that provides simple tools for doing phylogenetics.',
    author='Zach Sailer',
    author_email='zachsailer@gmail.com',
    packages=['phylogenetics'],
    install_requires=[
        'numpy',
    ],
    zip_safe=False)
