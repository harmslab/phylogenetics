import os
import subprocess

def runCdhit(homolog_list,redund_cutoff=0.99,tmp_file_suffix="oB_cdhit",
             keep_tmp=False):
    """
    Remove redundant homologs using cdhit.  After clustering with
    the redundancy cutoff, take the member of each cluster with the lowest
    rank.  Return a subset of homolog_list.
    """

    # Write out the fasta file with a unique name for each sequence that goes
    # >0, >1...>N.  Those numbers point to the index of the sequence in
    # homolog_list.

    # Don't do anything for empty list
    if len(homolog_list) == 0:
        print("Warning: empty list passed to cdhit!  Ignoring.")
        return homolog_list

    fasta_string = "".join([s.formatFasta(i) for i,s in enumerate(homolog_list)])
    f = open("%s.fasta" % tmp_file_suffix,'w')
    f.write(fasta_string)
    f.close()

    # Run cdhit
    cdhit_cmd = "cdhit -i %s.fasta -o %s_cdhit -c %.3f" % (tmp_file_suffix,
                                                           tmp_file_suffix,
                                                           redund_cutoff)
    args = shlex.split(cdhit_cmd)

    run = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdoutdata, stderrdata = run.communicate()
    if run.returncode != 0:
        print(stdoutdata)
        err = "cdhit failed!\n"
        raise BlastToolsError(err)

    # Now parse the output of cdhit and grab members of clusters with the
    # lowest rank
    f = open("%s_cdhit.clstr" % tmp_file_suffix,'r')

    out = []
    in_cluster = []
    line = f.readline()
    while line != "":

        # If we are starting a new cluster
        if line.startswith(">"):

            # ... and this is not the first cluster
            if in_cluster != []:

                # Take the member of in_cluster with the minimum rank
                ranks = [homolog_list[c].rank for c in in_cluster]
                best = in_cluster[ranks.index(min(ranks))]
                out.append(homolog_list[best])

            in_cluster = []

        # If this is not a blank line, record the seq_id in in_cluster
        elif line[0] in string.digits:
            seq_id = int(line.split(">")[1])
            in_cluster.append(seq_id)

        # Read the next line
        line = f.readline()

    # Grab the last cluster
    ranks = [homolog_list[c].rank for c in in_cluster]
    best = in_cluster[ranks.index(min(ranks))]
    out.append(homolog_list[best])

    f.close()

    # Delete temporary files
    if not keep_tmp:
        os.remove("%s.fasta" % tmp_file_suffix)
        os.remove("%s_cdhit" % tmp_file_suffix)
        os.remove("%s_cdhit.clstr" % tmp_file_suffix)
        os.remove("%s_cdhit.bak.clstr" % tmp_file_suffix)

    print("cdhit lowered redundancy @ %.3f, %i of %i kept" % (redund_cutoff,
                                                              len(out),
                                                              len(homolog_list)))

    return out



def cdhit_homologs(homologs):
    """ Run function that clusters homologs. """
    pass 