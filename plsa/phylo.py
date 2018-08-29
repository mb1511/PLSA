'''
Created on 31 Aug 2016

@author: Matt

@change: Add user options
@change: Overhauled module
'''
from __future__ import print_function

import subprocess as sp
import numpy as np
import os
import glob
from os.path import join, basename

import multi
import fasta
import f_con
import fmtdata
import muscle
from mega import distance, read_meg

BLAST_PATH = None # path to NCBI BLAST+ bin/ directory

@multi.map_list
def resample(sample_no=0, gene_list=[], wd='src'):
    '''
    Aligns and creates ML Tree from sampled genes
    
    Utilises mulitple processes
    '''
    tree = join(wd, 'sampling/tree%d.meg' % sample_no)
    aln_path = join(wd, 'sampling/aln%d.fasta' % sample_no)
    conc_path = join(wd, 'sampling/conc%d.fasta' % sample_no)
    # write concatenated gene list to file
    # run muslce and make tree
    with open(conc_path, 'w') as conc:
        conc.write(gene_list)
        # run alignment (using MUSCLE application)
        
    mur = muscle.run('-in', conc_path, '-maxiters', '2', path=aln_path, get_aln=False)
    # make NWK distance tree using MEGA7 commandline application
    if not mur:
        mtr = distance.run(aln_path, out=tree)
        #mtr = mt.max_likelihood(aln_path, bootstrap=False, boot_num=0, out=tree)
        # return tree path when made
        if not mtr:
            return tree
    return 'Sample %d Failed' % sample_no

def make_gene_list(seq_lists, size=10, repeats=1000, wd='src'):
    '''
    Generator Function to create samples
    
    seq_lists    - sequence of file paths to conserved gene lists
    size         - size of sample
    repeats      - number of samples to take
    '''
    # TODO: add check to see if samples dir is empty
    for rep in xrange(repeats):
        gene_list = ''
        start = True
        max_l = 0
        for faa in seq_lists:
            title = os.path.basename(faa)
            seqs = fasta.fasta_read(faa, generator=False)
            if not start:
                # should always be true - if not, something has gone
                # wrong somewhere
                assert len(seqs) == max_l
            else:
                max_l = len(seqs)
            if start:
                start = False
                # get some random indices
                indices = [i for i in np.random.choice(range(len(seqs)),
                                                       size,
                                                       replace=False)]
                for i in indices:
                    # write sample genes to file
                    seq = seqs[i]
                    #print(seq.name)
                    with open(join(wd, 'sampling/accs%d.txt' % rep), 'a') as s:
                        s.write(title + '\t'  + seq.name[1:seq.name.find(' ')] + '\n')
                            
            gene_list += '>' + title + '\n'
            gene_list += ''.join(seqs[i].seq for i in indices)
            gene_list += '\n'
        yield (rep, gene_list)

def make_single_genes(seq_lists, wd='src'):
    seqs = fasta.fasta_read(seq_lists[0], generator=False)
    for i in xrange(len(seqs)):
        start = True
        gene_list = ''
        for faa in seq_lists:
            title = os.path.basename(faa)
            seqs = fasta.fasta_read(faa, generator=False)
            if start:
                start = False
                with open(join(wd, 'sampling/accs%d.txt' % i), 'a') as s:
                    s.write(title + '\t'  + seqs[i].name[1:seqs[i].name.find(' ')] + '\n')
            gene_list += '>' + title + '\n'
            gene_list += seqs[i].seq
            gene_list += '\n'
        yield (i, gene_list)

def run(db_files=[], cov=65, size=10, make_dbs=True, fcon=True, repeats=500,
        wd='', default_genes='', num_threads=7, singles=False, ext='faa'):
    '''
    Perfroms phylogenetic analysis
    
    make blast data
       |
    find conserved
       |
    format results and save files
       |
    generate concatenated random samples
       |
    align sample file using MUSCLE
       |
    get pairwise distances using MEGA
    '''    
    
    while 1:
        cont = raw_input(
            'The following directories will be created:\n%s\n%s\n%s\nDo you wish to continue? (y/n): ' % (
                join(wd, 'blast_data'),
                join(wd, 'sampling'),
                join(wd, 'ref_genes'))
        if cont.lower() == 'y':
            break
        elif cont.lower() == 'n':
            return
        else:
            print('Please select "y" or "n"...')
    
    makeblastdb = join(BLAST_PATH, 'makeblastdb')
    
    if not os.path.isdir(join(wd, 'blast_data')):
        os.mkdir(join(wd, 'blast_data'))
        
    if not os.path.isdir(join(wd, 'sampling')):
        os.mkdir(join(wd, 'sampling'))
    
    # TODO: multiprocessing
    def make_db(faa, db_type='prot'):  
        sp.call([
            makeblastdb, '-in', faa,
            '-parse_seqids','-dbtype', db_type,
            '-out', join(wd, 'blast_data', basename(faa))])
    
    if make_dbs:
        # make all databases individually
        # clear data folder
        index = 0
        for index, f in enumerate(glob.glob(join(wd, 'blast_data/*.*'))):
            os.remove(f)
        print('%d files removed from blast_data/' % index)
        add.status('Folder cleared')
        # make blast data files
        for f in db_files:
            make_db(f)
            # copy over faa file to blast data directory
            with open(join(wd, 'blast_data', basename(f)), 'w') as s:
                with open(f) as r:
                    s.write(r.read())
        add.status('Databases constructed')
    
    if fcon:
        print('Finding conserved genes...')
        f_con.run_new(
            default_genes, join(wd, 'blast_data'),
            cov=cov, wd=wd, num_threads=num_threads)
     
    # TODO: add some status updates
    seq_lists = glob.glob(join(wd, 'ref_dbs/*.%s' % ext))
    
    if not singles:
        # create samples
        gene_lists = make_gene_list(seq_lists=seq_lists, size=size, repeats=repeats, wd=wd)
        dists = _fmtdata.Pairs()
        for tree in resample(arg_list=gene_lists, cores=num_threads, wd=wd):
            #for pair in _read_tree.get_distances(tree):
            for pair in _read_meg.get_all_dists(tree):
                try:
                    dists[pair[0], pair[1]]
                    dists[pair[0], pair[1]].append(pair[2])
                except KeyError:
                    dists[pair[0], pair[1]] = [pair[2]]
                        
        # save progress
        with open(join(wd, 'sampling/dists.txt'), 'w') as s:
            s.write('Samples: %d\n' % repeats)
            s.write('Strain 1\tStrain 2\tMean\tStandard Deviaiton\n')
            for g1, g2 in dists:
                x = dists[g1, g2]
                dist = str(np.mean(x))
                if len(x) > 1:
                    # catch div/0 error
                    sd = str(np.std(x))
                else:
                    # n = 0
                    sd = str(0.0)
                s.write(g1 + '\t' + g2 + '\t' + dist + '\t' + sd + '\n')
        
        # Return something
        with open(join(wd, 'sampling/dists.txt')) as d:
            out_text = d.read()
        return out_text
        
    else:
        gene_lists = make_single_genes(seq_lists=seq_lists, wd=wd)
        dists = _fmtdata.Pairs()
        for tree in resample(arg_list=gene_lists, cores=num_threads, wd=wd):
            for pair in _read_meg.get_all_dists(tree):
                id_ = '!REDACTED!'
                if id_ in pair[0] or id_ in pair[1]:
                    try:
                        dists[pair[0], pair[1]]
                        dists[pair[0], pair[1]].append(pair[2])
                    except KeyError:
                        dists[pair[0], pair[1]] = [pair[2]]
        seqs = fasta.fasta_read(join(wd, 'ref_dbs/!REDACTED!'))
        
        
        seq_ids = [seq.name[1:seq.name.find(' ')] for seq in seqs]
        out = '\t\t%s\n' % '\t'.join(str(x) for x in seq_ids)
        for key, item in dists.items():
            out += '%s\t%s\t%s\n' % (key[0], key[1], '\t'.join(str(x) for x in item))
        return out
            
            
    
   




