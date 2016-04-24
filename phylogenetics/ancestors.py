def unique_id(start, end=None):
    """ Returns a set of unique names """
    if end is None:
        end = int(start)
        start = 0
    return ["ZZ%08d" % i for i in range(start, end+1)]

class Ancestor:

    def __init__(self, id, Tree, **kwargs):

        ## Quality control name!
        self._Tree = Tree
        self.id = id

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
