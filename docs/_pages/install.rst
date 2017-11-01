Installation and setup
======================

Phylogenetics is still under heavy development and not ready for general use. The
best way to install currently is to clone the repo from Github and use ``pip`` to
install a development version.

::

    git clone https://github.com/Zsailer/phylogenetics
    cd phylogenetics
    pip install -e .


External tools
--------------

It's important that the following executables are callable from any location
for ``phylogenetics`` to work properly with external tools:

1. ``blastp`` and ``blastn``
2. ``cdhit``
3. ``msaprobs``
4. ``phyml``


NCBI-BLAST+ application and BLAST legacy code.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The NCBI-BLAST+ application is used to BLAST remote NCBI servers and local databases. Install this application by downloading the latest release of `NCBI-Blast`_ for your machine. Download and unpack.

.. _NCBI-Blast: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

For parts of this package, you may also need the NCBI legacy code (particularly for clustering with CDHIT under a threshold of 0.4). Download the `Legacy`_ code.

.. _Legacy: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/

Export the path of both these software packages.

Cdhit
^^^^^

CDhit is used to remove redundancy in a set of sequences for the reconstruction calculation. It clusters sequences that wouldn't add much new information to the maximum likelihood calculation.

TO install, clone the `CDHIT`_ repo and follow the simple instructions on this page (just `make` commands).

.. _CDHIT: https://github.com/weizhongli/cdhit

Export the path of to this installation in your ``bashrc`` file.

Sometimes the default executable to run is ``cd-hit`` rather than ``cdhit``. You may need to change this by changing the name of this executable to ``cdhit``.

Also, for ``cdhit`` runs with thresholds < 0.4, you must also export the path to ``psi-cd-hit.pl``. This is buried in ``cdhit``. This should just include an extra directory, ``psi-cd-hit``.

MSAProbs
^^^^^^^^
MSAProbs is used in this package to build a multiple sequence alignment. Download, unpack, and install the latest release of `MSAProbs`_.

.. _MSAProbs: http://sourceforge.net/projects/msaprobs/

Export the path to the MSAProbs executable to your ``bashrc``.

PhyML
^^^^^

Get the latest, stable release of `PhyML`_

.. _PhyML: https://github.com/stephaneguindon/phyml-downloads/releases

Download, unpack, and install Phyml. Export path to the ``src`` file in ``bashrc`` file.
