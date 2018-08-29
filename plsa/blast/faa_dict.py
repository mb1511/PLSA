'''
Created on 5 Sep 2016

@author: Matt
'''

from ... import fasta

class faa(dict):
    
    def __init__(self, faa_path):
        genes = fasta.fasta_read(faa_path)
        _dict = dict( [ (g.name[1:g.name.find('.')], g.seq) for g in genes   ] )
        super(faa, self).__init__(_dict)
    
    def __getitem__(self, key):
        try:
            ret = super(faa, self).__getitem__(key)
        except KeyError:
            key = key[:key.find('.')]
            ret = super(faa, self).__getitem__(key)
        return ret



