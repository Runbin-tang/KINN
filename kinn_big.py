#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 17:18:31 2021

@author: runbin
"""
from sys import exit
from itertools import product
from collections import Counter
import gc,random
import numpy as np
import pandas as pd

Protein=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

protein={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,
     'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}

dna={'A':0,'C':1,'G':2,'T':3}

def num2str(kmernum,k,base):
    kmer=list('A'*k)
    if base=='dna':
        bases={0:'A',1:'C',2:'G',3:'T'}
    elif base=='protein':
        bases={pj:pk for pk,pj in protein.items()}
        
    for i in range(k-1,-1,-1):
        if kmernum!=0:
            temp=kmernum%len(bases)
            kmer[i]=bases[temp]
            kmernum=kmernum//len(bases)
    kmer="".join(kmer)
    return kmer

def str2num(ks,base):
    m=0
    for i in ks:
        m=len(base)*m+base[i]
    return m

def gen_ks(k,base):
    ks={}
    #base={'A':0,'C':1,'G':2,'T':3} 
    if base=='protein':
        bases=protein
    elif base=='dna':
        bases=dna
    else:
        exit('generate keys error/')
    s=''.join(bases.keys())
    for i in product(s,repeat=k):
        m=0
        for j in i:
           m=len(bases)*m+bases[j]
        ks[''.join(i)]=m
    return ks

def relatv(j,fracm):
    nus={}
    for ji,jk in Counter(j).items():
        if jk not in nus.keys():
            nus[jk]=[]
        nus[jk].append(ji)
    if max(nus.keys())!=1:
        jmode=np.mean(nus[max(nus.keys())])
        #jnumber=len(nus[max(nus.keys())])*max(nus.keys())
    else:
        jmode=np.mean(j)
        #jnumber=1
    jj=np.array(j)
    jj=jj/fracm
    jm=jmode/fracm
    #jm=jm*jnumber
    #original
    return sum(jj*np.log2(jj*jm))    #sum(jj/jm-1) 
    # - original 
    #return -sum(jj*np.log2(jj*jm)) 

def kinn_pro(k,seqs):
    import re
    name=seqs[0]
    seq=seqs[1]
    ind={} # record the inner distance of ks-ks
    snk={} # refresh the position
    k2p={}# record the first and the end position of kmer
    for i in range(len(seq)-k+1):
        sks=seq[i:i+k]
        if bool(re.search(r'[^A,C,D,E,F,G,H,I,K,L,M,N,P,Q,S,T,V,W,Y]', sks)):
            continue
            #sks=re.sub(r'[^A,C,D,E,F,G,H,I,K,L,M,N,P,Q,S,T,V,W,Y]',
            #           random.choice(Protein),sks)
        if sks not in k2p.keys():
            k2p[sks]=[0,0]
            k2p[sks][0]=i
            snk[sks]=-1
        k2p[sks][1]=i
        #inn
        for j in k2p.keys():
            if snk[j]>=snk[sks] and snk[j]>-1:  
                tempk= str2num(j+sks,protein)
                if tempk not in ind.keys():
                    ind[tempk]=[]
                ind[tempk].append(i-snk[j])
        snk[sks]=i

    for ki in ind.keys():
        i=num2str(ki, 2*k, 'protein')
        if ind[ki]:
            fracm=k2p[i[k:]][-1]-k2p[i[:k]][0]
            ind[ki]=relatv(ind[ki], fracm)
        
    del k2p,snk
    gc.collect()
    return [name,ind]

def kinn_(k,seqs): 
    import re
    name=seqs[0]
    seq=seqs[1]
    #seq=re.sub(r'[^A,C,G,T]','',seqs[1])
    #ind=[0]*(4**(2*k)) # record the inner distance of ks-ks
    ind={}             # record the inner distance of ks-ks
    snk={} # refresh the position
    k2p={}# record the first and the end position of kmer
    for i in range(len(seq)-k+1):
        sks=seq[i:i+k]
        if bool(re.search(r'[^A,C,G,T]', sks)):
            continue
        if sks not in k2p.keys():
            k2p[sks]=[0,0]
            k2p[sks][0]=i
            snk[sks]=-1
        k2p[sks][1]=i
        #inn
        for j in k2p.keys():
            if snk[j]>=snk[sks] and snk[j]>-1:  
                tempk= str2num(j+sks,dna)
                
                if  tempk not in ind.keys():
                    ind[tempk]=[]
                ind[tempk].append(i-snk[j])
        snk[sks]=i
 
    for ki in ind.keys():
        i=num2str(ki, 2*k, 'dna')
        if ind[ki]:
            fracm=k2p[i[k:]][-1]-k2p[i[:k]][0]
            ind[ki]=relatv(ind[ki], fracm)
    del k2p,snk
    gc.collect()
    return [name,ind]


# seq='ACDCDADCCDC'
# seq1='ACTCTATCCTC'
# k=1

# # from Bio import SeqIO
# # seq=SeqIO.parse('data/peerdata/dna/hiv.fasta', 'fasta')
# # s=[[i.id.split('.')[0],str(i.seq)] for i in seq]
# # # seq=str(seq.seq)
# seq=['ref',seq]
# seq1=['ref',seq1]
# # ks=gen_ks(k,'protein')
# # seq=s[1]
# ss=kinn_pro(k, seq)
# ss1=kinn_(k, seq1)

# ss=[]
# for i in s:
#     ss.append(kinn_(k, ks,i))

# ks2={}
# for i in product('ACGD',repeat=k*2):
#     m=0
#     for j in i:
#        m=4*m+base[j]
#     ks2[''.join(i)]=m

