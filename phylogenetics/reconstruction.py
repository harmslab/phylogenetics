from phylogenetics.exttools.paml import CodeML

class Reconstructed(object):

    def __init__(self, HomologSet):
        """
        """
        self._HomologSet = HomologSet
        self._Tree = self._HomologSet.Tree

    def run(self, paml_job=CodeML):
        """ Run a reconstruction
        """
        seqfile = "asr-alignment.fasta"
        outfile = "asr-output"
        treefile = "asr-tree.nwk"

        self._Tree._DendroPyTree.write(path=treefile, schema="newick", suppress_internal_node_labels=True)
        self._HomologSet.Alignment.Write.fasta(fname=seqfile)

        self.paml_job = paml_job(
            seqfile=seqfile,
            outfile=outfile,
            treefile=treefile,
            fix_alpha=True,
            alpha=self._Tree.stats["Gamma shape parameter"],
        )

        self.paml_job.run()
