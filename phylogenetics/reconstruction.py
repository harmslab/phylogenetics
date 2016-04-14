
from phylogenetics.paml.paml import CodeML



class Reconstruction(object):

    def ___init__(self, Tree, paml_job=CodeML):
        """
        """
        self.Tree = Tree
        self.paml_job = CodeML()

    def run(self):
        """ Run a reconstruction
        """
        pass
