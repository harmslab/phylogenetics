
from phylogenetics.paml.paml import CodeML



class Ancestor:

    def __init__(self, Tree, data):

        ## Quality control name!
        self._Tree = Tree
        self._data = data

    @property
    def posteriors(self):
        return self._posteriors

    @property
    def mlsequence(self):
        return self._mlsequence


class AncestorSet(object):

    def __init__(self, Tree, Ancestors=None):
        """ Object that holds a set of ancestor object. """

        self._Tree = Tree

        if Ancestors is not None:
            for a in Ancestors:
                self.add(a)

    def add(self, Ancestor):
        """ Add Ancestor to set. """
        setattr(self, Ancestors.name, Ancestors)

    def rm(self, Ancestor):
        """ Remove the Ancestor from a set. """
        delattr(self, Ancestor.name)


class Reconstructed(object):

    def ___init__(self, Tree, paml_job=CodeML):
        """
        """
        self.Tree = Tree
        self.paml_job = CodeML()

    def run(self):
        """ Run a reconstruction
        """
        pass
