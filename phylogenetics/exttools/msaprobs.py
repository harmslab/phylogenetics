# Multiple sequence alignment using MSAProbs
from __future__ import absolute_import

import os, shlex
import subprocess

def run(fasta_fname="alignment", cores=2, rm_tmp=True):
    """ Wrapper for MSAProbs.

        Runs a multiple sequence alignment for sequences in a named fasta file.

        Args:
        ----
        fasta_fname: str
            filename for temporary files
        rm_tmp: bool
            If true, removes the temporary files that it creates.

        Returns:
        -------
        oname : str
            Output fasta filename
    """

    # Create a temporary fasta file from homologs as input to CDHIT.
    fname = "%s.fasta" % fasta_fname
    oname = "%s-finished.fasta" % fasta_fname

    # Build msaprobs command
    msaprobs_cmd = "msaprobs %s.fasta -o %s-finished.fasta -num_threads %s" % (fasta_fname,
                                                           fasta_fname, str(cores))
    # Format the command into list.
    args = shlex.split(msaprobs_cmd)

    # Run msaprobs using args.
    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdoutdata, stderrdata = run.communicate()

    # Check if alignment worked correctly.
    if run.returncode != 0:
        print(stdoutdata)
        err = "msaprobs failed!\n"
        raise Exception(err)

    # Remove alignment file if desired
    if rm_tmp:
        os.remove("%s.fasta" % fasta_fname)

    return oname
