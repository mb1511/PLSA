'''
Created on 20 Nov 2017

@author: mb1511
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary:
'''
from __future__ import print_function

import subprocess as sp

import read_meg

MEGA = None
DIST = 'distance_estimation_pairwise_protein.mao'
    
def run(aln, out):
    
    p = sp.Popen(
        [MEGA, '-a', DIST, '-d', aln, '-o', out, '-n'],
        stdout=sp.PIPE, stderr=sp.PIPE)
    
    p.communicate()
    return p.returncode