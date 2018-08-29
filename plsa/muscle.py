'''
Created on 31 Aug 2016

@author: Matt
'''
from __future__ import print_function

import subprocess as sp

def run(*args, **kwargs):
    args = list(args)
    path = kwargs.pop('path', 'src/muscle/temp_aln.fasta')
    get_aln = kwargs.pop('get_aln', True)
    
    if '-out' not in args:
        args.extend( ['-out', path])
    else:
        x = args.index('-out')
        path = args[x+1]
    
    print('Running MUSCLE on %s' % args[1])
    m = sp.Popen([muscle] + args, stdout=sp.PIPE, stderr=sp.PIPE)
    o, e = m.communicate()
    rc = m.returncode
    
    if not rc:
        if get_aln:
            with open(path) as s:
                text = s.read()
            return text
        else:
            print('MUSCLE Alignment of %s complete' % (args[1]))
            return 0
    else:
        print('MUSLCE Failed\nSTDOUT:\n%s\nSTDERR:\n%s' % (o, e))
        return 1