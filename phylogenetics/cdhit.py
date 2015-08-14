import os, shlex, string 
import subprocess

from phylogenetics.base import rank_homologs

def run_cdhit(homolog_set,redund_cutoff=0.99,tmp_file_suffix="oB_cdhit",
             keep_tmp=False, positive=(), negative=("putative","hypothetical",
             "unnamed", "possible", "predicted", "unknown", "uncharacterized",
             "mutant", "isoform")):
    """
    Remove redundant homologs using cdhit.  After clustering with
    the redundancy cutoff, take the member of each cluster with the lowest
    rank.  Return the subset of homologlist.
    """

    # Write out the fasta file with a unique name for each sequence that goes
    # >0, >1...>N.  Those numbers point to the index of the sequence in
    # homolog_list.

    # Don't do anything for empty list
    if len(homolog_set.homologs) == 0:
        print("Warning: empty list passed to cdhit!  Ignoring.")
        return homolog_set

    # Ranks homologs.
    rank_homologs(homolog_set, positive=positive, negative=negative)

    # Create a temporary fasta file from homologs as input to CDHIT.
    fname = "%s.fasta" % tmp_file_suffix
    homolog_set.write(fname, format="fasta", tags=["id"])

    # Build cdhit command
    cdhit_cmd = "cdhit -i %s.fasta -o %s_cdhit -c %.3f" % (tmp_file_suffix,
                                                           tmp_file_suffix,
                                                           redund_cutoff)
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

    # Now parse the output of cdhit and grab members of clusters with the
    # lowest rank
    f = open("%s_cdhit.clstr" % tmp_file_suffix,'r')

    # Get a id-to-rank mapping dict from homolog_set
    id_rank = homolog_set.get_map("id", "rank")

    subset_ids = []
    in_cluster = []
    line = f.readline()
    while line != "":

        # If we are starting a new cluster
        if line.startswith(">"):

            # ... and this is not the first cluster
            if in_cluster != []:

                # Take the member of in_cluster with the minimum rank
                ranks = [id_rank[homolog_id] for homolog_id in in_cluster]
                best_id = in_cluster[ranks.index(min(ranks))]
                subset_ids.append(best_id)

            in_cluster = []

        # If this is not a blank line, record the seq_id in in_cluster
        elif line[0] in string.digits:
            seq_id = line.split(">")[1][:10]
            in_cluster.append(seq_id)

        # Read the next line
        line = f.readline()

    # Grab the last cluster
    ranks = [id_rank[homolog_id] for homolog_id in in_cluster]
    best_id = in_cluster[ranks.index(min(ranks))]
    subset_ids.append(best_id)
    f.close()

    # Delete temporary files
    if not keep_tmp:
        try:
            os.remove("%s.fasta" % tmp_file_suffix)
            os.remove("%s_cdhit" % tmp_file_suffix)
            os.remove("%s_cdhit.clstr" % tmp_file_suffix)
            os.remove("%s_cdhit.bak.clstr" % tmp_file_suffix)
        except:
            pass

    homolog_subset = homolog_set.subset_homologs(subset_ids)

    return homolog_subset
