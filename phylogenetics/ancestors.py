from collections import OrderedDict

from phylogenetics.dataio import ancestorio, ancestorsetio

def unique_id(start, end=None):
    """ Returns a set of unique names """
    if end is None:
        end = int(start)
        start = 0
    return ["ZZ%08d" % i for i in range(start, end+1)]

class Ancestor:

    def __init__(self, id, Tree):
        """
            Object that holds an Ancestors' node data
        """
        self._Tree = Tree
        self.id = id

    @property
    def ml(self):
        """ Get maximum likelihood data. """
        ml_data = []
        for site_i, site_data in self.sites.items():
            # Iterate through residues to find the most probable
            ml_site = ("A", 0.0)
            for residue, prob in site_data.items():
                if prob > ml_site[1]:
                    ml_site = (residue, prob)

            # Add ml_site to ml_data
            ml_data.append(ml_site)
        # Return as ordered dict
        return ml_data

    @property
    def posteriors(self):
        return [site[1] for site in self.ml]

    @property
    def mlsequence(self):
        """ Get the maximum likelihood sequence for the ancestor. """
        # Iterate through sites data and find highest ML
        return [site[0] for site in self.ml]

    def site(self, site_n):
        """ Return the posterior probabilities of each residue for a given site."""
        return self.sites[site_n]

class AncestorSet(object):

    def __init__(self, Tree, ancestors=[]):
        """ Object that holds a set of ancestor object. """

        self._Tree = Tree

        if ancestors is not None:
            for a in ancestors:
                self.add(a)

        self.Read = ancestorsetio.Read(self)
        self.Write = ancestorsetio.Write(self)

    def add(self, Ancestor):
        """ Add Ancestor to set. """
        setattr(self, Ancestor.id, Ancestor)

    def rm(self, Ancestor):
        """ Remove the Ancestor from a set. """
        delattr(self, Ancestor.id)

    @property
    def ancestors(self):
        return self._ancestors
