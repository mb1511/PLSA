'''
Created on 20 Nov 2017

@author: mb1511
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary:
'''

import re
import fmtdata

def read_taxa(f):
    s_flag = False
    with open(f) as s:
        for line in s:
            if line.startswith('['):
                s_flag = True
                yield line[line.find('#')+1:line.find('\n')]
            else:
                if s_flag:
                    break

def get_dists(f):
    s_flag = False
    d_flag = 0
    with open(f) as s:
        for line in s:
            if not d_flag:
                if line.startswith('['):
                    s_flag = True
                else:
                    if s_flag:
                        d_flag = 1
            elif line.startswith('['):
                if d_flag == 1:
                    d_flag += 1
                    continue
                elif d_flag == 2:
                    yield [0.0]
                    d_flag += 1
                    continue
                else:
                    vals = map(float, re.findall(r'[+-]?[0-9.]+', line[line.find(']')+1:]))
                    if vals:
                        yield vals + [0.0]
                    

def get_all_dists(f):
    taxa = [t for t in read_taxa(f)]
    for i, vals in enumerate(get_dists(f)):
        for j, val in enumerate(vals):
            yield (taxa[i], taxa[j], val)

def gen_pairs(f):
    p = _fmtdata.Pairs()
    for x,y,z in get_all_dists(f):
        p[x, y] = z
    return p

def get_avg(f):
    tot = 0
    n = 0
    for x,y,z in get_all_dists(f):
        tot += z
        n += 1.
    return tot/n
        

if __name__ == '__main__':
    import glob
    from os.path import basename, splitext
    accs = {}
    with open('../../ref_dbs/GCA_000007325.1_ASM732v1_protein.faa') as s:
        for line in s:
            if line.startswith('>'):
                accs[line[1:9]] = []
    
    
    for sample in glob.glob('../../sampling/*.meg'):
        print sample
        avg = get_avg(sample)
        snam = splitext(basename(sample))[0]
        n = snam[4:]
        ac_list = []
        with open('../../sampling/accs%s.txt' % n) as r:
            for line in r:
                if line.startswith('G'):
                    ac_list.append(line.split('\t')[1][:-1])
        for acc in ac_list:
            accs[acc].append(avg)
    for k, i in accs.items():
        print k, float(sum(i))/float(len(i))
    
    
    
    
    
