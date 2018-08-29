'''
Created on 12 Sep 2017

@author: mb1511
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary: Read NWK trees
'''

from ete3 import Tree

def _get_leaves(t):
    for node in t.traverse('postorder'):
        if node.is_leaf():
            yield node

# leaf matirx
def get_distances(tree_path):
    t = Tree(tree_path)
    for leaf_1 in _get_leaves(t):
        for leaf_2 in _get_leaves(t):
            yield (leaf_1.name, leaf_2.name, leaf_1.get_distance(leaf_2))