import numpy as np
import json
import pickle


class Homolog(object):

    def __init__(self, unique_id, sequence, **kwargs):
        
        # Must set a unique ID and sequence
        self.id = unique_id
        self.sequence = sequence
        
        # Set the other attributes of 
        print(kwargs)
        for key, value in kwargs.items():
            setattr(self, key, value)
    
    def mapping(self, key1, key2):
        """ Return mapping between two attributes in homolog """
        return dict(getattr(self, key1), getattr(self, key2))
    
    # -----------------------------------
    # Output formats
    # -----------------------------------
    
    def fasta(self,tags=None):
        """ Return fasta formatted string with named tags (in order given).
        
            If no tags are given, prints id with sequence.
        """
        f = ">"
        if tags is not None:
            # Remove sequence if its an elements
            if "sequence" in tags:
                tags.remove("sequence")
            f += "|".join([str(getattr(self, t)) for t in tags])
        else:
            f += self.id
        # New line with full sequence string
        f += "\n" + self.sequence
        return f
        
    def json(self, **kwargs): 
        """ Return json formatted string. """
        return json.dumps(self.__dict__)
    
    def pickle(self, **kwargs):
        """ Returns pickle string. """
        return pickle.dumps(self)
        
    def write(self, filename, format="fasta", tags=None):
        """ Write to file with given format. 
        
            Default is fasta.
        """
        format_func = getattr(self, format)
        f = open(filename, "w")
        f.write(format_func(tags))
        f.close()
        
        
class HomologSet(object):
    
    def __init__(self, homolog_set=[]):
        """ Construct a set of homolog objects. """
        self._homologs = homolog_set
        
    @property
    def homologs(self):
        """ Get homolog set. """
        return self._homologs
    
    def mapping(self, attr1, attr2):
        """ Return mapping between two attributes in homolog set."""
        m = dict()
        for h in self._homologs:
            m[getattr(h, attr1)] = getattr(h, attr2)
        return m
    
    def add_homolog(self, homolog):
        """ Append a homolog object to the set."""
        if isinstance(homolog, Homolog): 
            self._homologs.append(homolog)
        else:
            raise Exception("homolog must be an instance of Homolog class.")
                    
    def fasta(self, tags=None):
        """ Return string in fasta format for the set."""
        f = ""
        for h in self._homologs:
            f += h.fasta(tags) + "\n"
        return f
    
    def json(self, **kwargs):
        """ Return json string of homolog set."""
        obj = list()
        for h in self._homologs:
            obj.append(h.__dict__)
        return json.dumps(obj)
    
    def pickle(self, **kwargs):
        """ Return pickle string of homolog set. """
        return pickle.dumps(self)
    
    def write(self, filename, format="fasta", tags=None):
        """ Write to file with given format. 
        
            Default is fasta.
        """
        format_func = getattr(self, format)
        f = open(filename, "w")
        f.write(format_func(tags))
        f.close()