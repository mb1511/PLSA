'''
Created on 10 Aug 2016

@author: Matt
'''

import os

import multi
import fasta
from blast import local_blast, faa_dict
    
def compare(database, file_path=None, from_text=None, from_list=None, cov=65, remove=[], **kw):
    removed = 0      
    
    # implemented for multi blast
    prot = _faa_dict.faa(database)
    if not from_list:
        g = fasta.fasta_read(file_path=file_path, from_text=from_text, generator=False)
    else:
        g = from_list
    print len(g)
    ref = []
    # create blast query
    with open('src/tools/blast/query.fasta', 'w') as q:
        q.write('\n'.join(gene.fasta for gene in g))
    # run (multi) blast search
    aln = local_blast.run(db_path=database, mr=1, make_db=False,
                          r_a=True, outfmt='details', join_hsps=True, **kw)
    
    for index, gene in enumerate(g):
        try:
            # grab top alignment (position 0) % similarity
            a = aln[index][0]
            sim = a.psim
        except IndexError:
            # if no record, move onto the next gene
            for j in xrange(len(remove)):
                del remove[j][index - removed]
                removed += 1
            continue
    
        if sim < cov:
            # if psim lower than cut off value, move onto next gene
            for j in xrange(len(remove)):
                del remove[j][index - removed]
                removed += 1
            continue
        else:
            print index, gene.name[1:gene.name.find(' ')], a.h_acc, sim
            # matched lists
            ref.append(fasta.sequence(a.h_def, prot[a.h_acc]))
            # keep track of genes to keep from previous queries
            #indices.append(index)
    
    remove.append(ref)
    return remove

@multi.map_list
def comps(thread, database, cov=65, ref=None):
    ret = []
    prot = _faa_dict.faa(database)
    query = 'src/tools/blast/query%d.fasta' % thread
    # create blast query
    with open(query, 'w') as q:
        q.write('\n'.join(gene.fasta for gene in ref))
    # run (multi) blast search
    aln = local_blast.run(
        db_path=database, query=query, mr=1, make_db=False,
        r_a=True, outfmt='details', join_hsps=True)
    for index, gene in enumerate(ref):
        try:
            # grab top alignment (position 0) % similarity
            a = aln[index][0]
            sim = a.psim
        except IndexError:
            # if no record, move onto the next gene
            ret.append((gene.index, None))
            continue
    
        if sim < cov:
            # if psim lower than cut off value, move onto next gene
            ret.append((gene.index, None))
            continue
        else:
            ret.append((gene.index, fasta.sequence(a.h_def, prot[a.h_acc], index=gene.index)))
            #print index, gene.name[1:gene.name.find(' ')], a.h_acc, sim
    return ret

def run_new(start='', d='blast_data', cov=65, wd='', **kw):
    
    num_threads = kw.pop('num_threads', 7)
    dbs = [str(d + '/' + f[:-4]) for f in os.listdir(d) if f.endswith('.phr')]   
    ref_genes = fasta.fasta_read(start, generator=False)
    matches = {}
    
    if not os.path.isdir(os.path.join(wd, 'ref_dbs')):
        os.mkdir(os.path.join(wd, 'ref_dbs'))
    
    # segment dbs into blocks (size = num_threads) - blocks don't have to be full
    # process blocks
    # refine
    # repeat until all blocks completed
    blocks = []
    for i in xrange(0, len(dbs), num_threads):
        blocks.append(dbs[i:i + num_threads])
    
    for block in blocks:
        args = zip(range(len(block)), block)
        to_remove = []
        print len(ref_genes)
        for i, res in enumerate(
            comps(
                arg_list=args,
                cov=cov,
                ref=ref_genes,
                cores=num_threads)):
            for r in res:
                if r[1] is None:
                    to_remove.append(r[0])
                else:
                    try:
                        matches[r[0]][block[i]] = r[1]
                    except KeyError:
                        matches[r[0]] = {block[i]: r[1]}
        # remove 0 score proteins from master list
        for r in set(to_remove):
            print matches.pop(r)[dbs[5]].index, 'Removed'
        ref_genes = [matches[g][dbs[5]] for g in matches]
    
    for d in dbs:
        with open(os.path.join(wd, 'ref_dbs', os.path.basename(d)), 'w') as s:
            for g in matches:
                s.write(matches[g][d].fasta)
                s.write('\n')
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    