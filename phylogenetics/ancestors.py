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
        self._gaps = []

    def _set_gaps(self, indices):
        """ Indices that should be gaps. """
        self._gaps = indices
        setattr(self, "gaps", self._return_gaps)

    @property
    def _return_gaps(self):
        return self._gaps

    @property
    def ml(self):
        """ Get maximum likelihood data. """
        ml_data = []

        for site_i, site_data in self.sites.items():
            # Iterate through residues to find the most probable
            ml_site = ("A", 0.0)

            # Check if a gap should be placed
            if site_i in self._gaps:
                ml_site = ("-", "-")

            # If no gap, find ML for that site
            else:
                for residue, prob in site_data.items():
                    if prob > ml_site[1]:
                        ml_site = (residue, prob)

            # Add ml_site to ml_data
            ml_data.append(ml_site)
        # Return as ordered dict
        return ml_data

    @property
    def site_pp(self):
        """ Return site-by-site posterior probability. """
        return [site[1] for site in self.ml]

    @property
    def mlsequence(self):
        """ Get the maximum likelihood sequence for the ancestor. """
        # Iterate through sites data and find highest ML
        return [site[0] for site in self.ml]

    @property
    def posterior(self):
        """ Return the average posterior probability of the ancestor."""
        vals = [val for val in self.site_pp if val != "-"]
        return sum(vals)/len(vals)

    def site(self, site_n):
        """ Return the posterior probabilities of each residue for a given site."""
        return self.sites[site_n]


class AncestorSet(object):

    def __init__(self, Tree, ancestors=[]):
        """ Object that holds a set of ancestor object. """

        self._Tree = Tree
        self._ancestors = {}
        if ancestors is not None:
            for a in ancestors:
                self.add(a)

        self.Read = ancestorsetio.Read(self)
        self.Write = ancestorsetio.Write(self)

    def _nodes_to_ancestor(self):
        """ Point nodes in tree to ancestor object and vice versa."""
        # Iterate through nodes in tree.
        for node in self._Tree._DendroPyTree.nodes():

            # Find nodes that represent tips of the tree.
            if node.taxon is None:
                # Get the id of the homolog
                id = node.label

                # Get Homolog object for this node
                Ancestor = getattr(self, id)

                # Bind Homolog to node
                node.Ancestor = Ancestor

                # Bind node object to Homolog
                setattr(Ancestor, "node", node)


    def add(self, Ancestor):
        """ Add Ancestor to set. """
        setattr(self, Ancestor.id, Ancestor)
        self._ancestors[Ancestor.id] = Ancestor

    def rm(self, id):
        """ Remove the Ancestor from a set. """
        delattr(self, id)
        del self._ancestors[id]

    @property
    def ancestors(self):
        return self._ancestors
