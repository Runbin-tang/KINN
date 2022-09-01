#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 11:02:34 2022

@author: runbin
"""

'''
need to compute the distance of value based key in dict
'''

import pandas as pd
import numpy as np
import tqdm
# from multiprocessing import Pool as pl
# from scipy.spatial.distance import pdist,squareform
# from functools import partial


def distance_dict(d):
    #pla=pl()
    
    dis=np.zeros((len(d),len(d)))
    for i in tqdm.tqdm(range(1,len(d))):
        #st=di(d,i)
        df = pd.DataFrame([d[i]]).T
        df = df.reset_index().rename(columns={'index':'key'})
        tp=np.array(list(d[i].values()))
        #js=[];
        cl=np.sqrt(sum(tp**2))
        df[0]=df[0]/cl
        for j in range(i):   
            tp=np.array(list(d[j].values()))
            df[j+1] = df['key'].apply(lambda x: d[j].get(x,0))
            df[j+1] = df[j+1]/np.sqrt(np.sum(tp**2))
            #js.append(cl*np.sqrt(np.sum(tp**2)))
        ns=df.loc[:, df.columns != 'key'].values
        ns=np.dot(ns[:,0],ns[:,1:])
        st=1-ns
        for j in range(i):
            dis[j][i]=dis[i][j]=st[j]
    return dis            
    