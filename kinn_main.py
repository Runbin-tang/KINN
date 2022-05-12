#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 08:55:53 2021

@author: runbin
"""
from scipy.spatial.distance import pdist,squareform
from multiprocessing import Pool as pl
from functools import partial
from collections import Counter
from Bio import SeqIO

import numpy as np
import pandas as pd
import math, argparse, gc

Protein=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

protein={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,
     'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}


def write_megaa(path,name,dist):
    
    with open(path,"w") as f:   
        f.write('#mega\n!TITLE;\n')
        for namei in name: 
            f.write('#'+namei+'\n')
        l=dist.shape[0]
        for i in range(l):
            for j in range(i):
                f.write(str('{:.10f}'.format(dist[i][j]))+'\t')
            f.write('\n')
    f.close()
    
def kinn(opt, seq, name):
    import kinn_com,tqdm
    tempcheck=1
    print('[Step 1]  obtain the  kinn of sequence (k={})'.format(opt.k))

    pa=pl()
    if opt.seqtype=='dna':
        
        if opt.k<7:
            ks=kinn_com.gen_ks(opt.k,opt.seqtype)
            pas=partial(kinn_com.kinn_, opt.k, ks)
            inn_f = list(tqdm.tqdm(pa.imap_unordered(pas,seq), 
                           desc='  -Progress', 
                           total=len(seq),ncols=100))
        else:
            import kinn_big
            
            pas=partial(kinn_big.kinn_, opt.k)
            inn_f = list(tqdm.tqdm(pa.imap_unordered(pas,seq), 
                           desc='  -Progress', 
                           total=len(seq),ncols=100))
            tempcheck=0
    elif opt.seqtype=='protein':
        
        if opt.k<4:
            ks=kinn_com.gen_ks(opt.k,opt.seqtype)
            pas=partial(kinn_com.kinn_pro, opt.k, ks)
            inn_f = list(tqdm.tqdm(pa.imap_unordered(pas,seq), 
                           desc='  -Progress', 
                           total=len(seq),ncols=100))
        else:
            import kinn_big
            pas=partial(kinn_big.kinn_pro, opt.k)
            inn_f = list(tqdm.tqdm(pa.imap_unordered(pas,seq), 
                           desc='  -Progress', 
                           total=len(seq),ncols=100))
            tempcheck=0
    else:
        exit('the sequence type is error/')
    
    print('[Step 2]  Compute distance matricx')
    names=[];ns=[]
    for i in inn_f:
        names.append(i[0])
        ns.append(i[1])
    del inn_f
    if tempcheck:
        distA=pdist(ns,metric=opt.metric)
        distB = squareform(distA) 
    else:
        import dic_distance
        distB = dic_distance.distance_dict(ns)
    
    print('[Step 3]  wrirte distance matrix into file.')
    if name:
        names=[name[i] for i in names]

    savefile=str(opt.k)+opt.seqtype +'_'+splitn(opt.seqs)
    path=os.path.join(opt.savefolder,savefile)
    write_megaa(path,names,distB)
    del distA, distB
    gc.collect()
    #return inn_f,distB

def splitn(s):
    if s.split('/')[-1]=='':
        return s.split('/')[-2].split('.fasta')[0]+'.meg'
    else:
        return s.split('/')[-1].split('.fasta')[0]+'.meg'
       
def optil(s):
    ave=0
    for i in range(len(s)):
        ave+=len(s[i][1])
    return math.log10(ave/len(s))-0.5

def parameters():
    
    parser = argparse.ArgumentParser()
        
    parser.add_argument('seqtype', type=str, choices=['dna','protein'], 
                        help=' the sequence stype')
    parser.add_argument('--k',  type=int, help='kmer')
    parser.add_argument('--seqs', default='../data_sample/sixteen_srna.fasta', 
                        type=str, help='the seqs file (all sequence in one fasta file)')
    parser.add_argument('--metric', default='cosine', type=str, 
                        help='metric of distance')
    parser.add_argument('--seqformat', default='fasta', type=str, 
                        help='the data stype')
    parser.add_argument('--medatacsv', type=str,
                        help='the sequence information')
    parser.add_argument('--savefolder', default='./distance', type=str, 
                        help='position of save distance file')
    
    opt = parser.parse_args()

    return opt

if __name__=='__main__':

    import os
    opt=parameters()
    data=os.path.join(opt.seqs)
    
    print('\n','='*8,'START ','='*8)
    print('[Step 0]  upload data')
    
    if os.path.isfile(data):
        s=[];recordname={}
  
        seq=SeqIO.parse(data, opt.seqformat)  
        s=[[i.id.split('.')[0],str(i.seq)] for i in seq]
        
        if opt.medatacsv:
            cs =pd.read_csv(opt.medatacsv,sep='\t',header=None)
            recordname=cs.set_index(cs[0])[1].to_dict()
            recordname={i.split('.')[0]:j for i,j in recordname.items()}
        else:
            recordname={s[i][0]:s[i][0] for i in range(len(s))}
            
    elif os.path.isdir(data):
        s=[];recordname={}
        fs=list(os.walk(data))[0][2]
        for i in fs:
            seq=SeqIO.read(os.path.join(data,i),opt.seqformat)
            s.append([seq.id, str(seq.seq)])
            recordname[seq.id]=i.split('.')[0]
    else:
        exit('data upload error, check your input !')
    
    if not opt.k:
        if opt.seqtype=='dna':
            opt.k=round(optil(s)/math.log10(4))
        else:
            opt.k=round(optil(s)/math.log10(20))
            
    if not os.path.exists(opt.savefolder):
        os.makedirs(opt.savefolder)
    print(' -sequence type is {}'.format(opt.seqtype))
    print(' -total {} sequences'.format(len(s)))
    print(' -distance metric is {}'.format(opt.metric))
    print(' -save distance mega file pathe: {}'.format(
                            os.path.join(opt.savefolder,str(opt.k)
                                   +opt.seqtype+'_'+splitn(opt.seqs))))

    ss=kinn(opt,s,recordname)
        
    print('='*8,'ENDING','='*8,'\n')
