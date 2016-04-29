import os, shlex, string
import subprocess

def run(fname_prefix,
    redund_cutoff=0.99,
    tmp_file_suffix="oB_cdhit",
    word_size=5,
    cores=1,
    keep_tmp=False,
    accession=(),
    positive=(),
    negative=("putative","hypothetical", "unnamed", "possible", "predicted",
                "unknown", "uncharacterized","mutant", "isoform")):
    """
        Remove redundant homologs using cdhit.  After clustering with
        the redundancy cutoff, take the member of each cluster with the lowest
        rank.  Return the subset of homologlist.
    """
    # Build cdhit command
    # If the threshold is below 0.4, must use psi-cd-hit command. This requires
    # blast legacy to be installed and the psi-cd-hit directory be exported to PATH
    if redund_cutoff < 0.4:
        cdhit_cmd = "psi-cd-hit.pl -i %s.fasta -o %s_cdhit -c %.3f -core %d" % (fname_prefix,
                                                               fname_prefix,
                                                               redund_cutoff,
                                                               cores)
    # Above the 0.4 threshold, we can use the typical cdhit
    else:
        cdhit_cmd = "cdhit -i %s.fasta -o %s_cdhit -c %.3f -n %d -T %d" % (fname_prefix,
                                                            fname_prefix,
                                                            redund_cutoff,
                                                            word_size,
                                                            cores)
    # Format the command into list.
    args = shlex.split(cdhit_cmd)

    # Run CDhit using args.
    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdoutdata, stderrdata = run.communicate()

    # Check if CDHIT worked correctly.
    if run.returncode != 0:
        print(stdoutdata)
        err = "cdhit failed!\n"
        raise BlastToolsError(err)
