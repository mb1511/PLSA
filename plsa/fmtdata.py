'''
Created on 29 Sep 2017

@author: mb1511
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary: organise distance data into different file formats
'''

class Pairs(dict):
    '''
    Class for creating a dict with unbiased tuple key order
    '''
    def __init__(self, *args):
        dict.__init__(self, *args)
    def __getitem__(self, key):
        key = tuple(sorted(key))
        return dict.__getitem__(self, key)
    def __setitem__(self, key, val):
        key = tuple(sorted(key))
        return dict.__setitem__(self, key, val)
    
class Nexus:
    def __init__(self, pairs):
        self.taxa = []
        self.pairs = pairs
        for n1, n2 in pairs.keys():
            if n1 not in self.taxa:
                self.taxa.append(n1)
            if n2 not in self.taxa:
                self.taxa.append(n2)
        self.taxa_labels = ["[%d]\t'%s'" % (i+1, n) for i, n in enumerate(self.taxa)]
        self.n_taxa = len(self.taxa)
        out= ''
        out += '#nexus\n'
        out += 'BEGIN Taxa;\n'
        out += 'DIMENSIONS ntax=%d;\n' % self.n_taxa
        out += 'TAXLABELS\n'
        out += '\n'.join(self.taxa_labels)
        out += '\n;END; [Taxa]\n'
        out += 'BEGIN Distances;\n'
        out += 'DIMENSIONS ntax=%d;\n' % self.n_taxa
        out += 'FORMAT labels=left diagonal triangle=lower;\n'
        out += 'MATRIX\n'
        for i, n in enumerate(self.taxa):
            out += self.taxa_labels[i] + '\t'
            out += '\t'.join(str(pairs[n, m]) if j <= i else '' for j, m in enumerate(self.taxa)) + '\n'
        out += ';\nEND; [Distances]'
        self.out = out
    
    def write(self, path):
        with open(path, 'w') as s:
            s.write(self.out)
        
        
        