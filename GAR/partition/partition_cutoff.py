# -*- coding: utf-8 -*-
"""
Partition the data for Peru using three cutoffs
rlafarguette@imf.org
Time-stamp: "2018-01-08 19:38:27 rlafarguette"
"""

###############################################################################
#%% Modules
###############################################################################
## Core modules
import os, sys, importlib
import pandas as pd          
import numpy as np

## Functions loading
from .partition import Partition
from datetime import datetime as date                 ## Dates

## Suppress warnings
import warnings
warnings.filterwarnings("ignore")


###############################################################################
# Fuction to do partion for one cutoff time. Return partion and loading
###############################################################################
def add_id(df, id_dict):
    """ Add identifiers variables to a pandas frame """
    variables_id = sorted(list(id_dict.keys()))
    for v, var in enumerate(variables_id):
        df.insert(v, var, id_dict[var])
    return(df)

def p_cutoff(dall,groups_dict,cutoff,bench,method, benchcutoff, PLStarget, saveim=False):


    df = dall.loc[dall.date >= cutoff].copy()
    
    #print(df)
    partition_dict = groups_dict
    variables=[]
    label_dict={}
    

    for key, values in partition_dict.items():
        variables.extend(values)
        for e in values:
            label_dict[e]=e
    

    #country=df['country'].values[0]
    c_id_dict = {'cutoff' : cutoff.strftime("%Y-%m-%d"), 
                 'variables': repr(variables)}
                #'variables': repr(gv.indvars_dict[country])}
    try:
        df.loc[:,'var_benchmark'] = df[bench]

        dfbch = df.dropna(subset=['var_benchmark'], axis=0, how='any')[['var_benchmark']].copy()
        ## Define the benchmark for the partition
        threshold = dfbch['var_benchmark'].quantile(benchcutoff) # Lower growth regime
        df.loc[:,'benchmark'] = (dfbch['var_benchmark'] < threshold)
        df['benchmark']=df['benchmark'].fillna(method='ffill')

        #print (df['benchmark'])
        #################### Run the partitionning on a subset
        #print(method)
        if method=='LDA':
            p = Partition(df, partition_dict,
                          reduction='LDA', benchmark='benchmark', PLStarget=PLStarget)
        elif method=='PCA':
            p = Partition(df, partition_dict,
                          reduction='PCA', benchmark=None , PLStarget=PLStarget)
        elif method=='PLS':
            p = Partition(df, partition_dict,
                          reduction='PLS', benchmark=None , PLStarget=PLStarget)
            

        dp = p.partition # Run the partition on the full frame

        
        dp = add_id(dp, c_id_dict)
        dp.loc[:,'date'] = dp.index
        dp.index.name=None

        ## Merge with the original data
#        id_cols = ['date']
#        old_vars = [x for x in df if x not in dp.columns] + id_cols
#        dp_all = pd.merge(dp, df[old_vars], on=id_cols, how='left')
#       dp_all = dp_all.set_index(dp_all.date, drop=False)
#        if saveim:
#            print(dp.values)
#            print('HERE')
#            print(df.values)
#            print(dp.columns)
#            print(dp_all.columns)
#            print(dp_all.values)
#            df.to_csv('df.csv')
#            dp.to_csv('dp.csv')
        ## Loading from the partitioning
        dl = p.loading; dl = add_id(dl, c_id_dict)
        dl.insert(0, 'variable_o', dl.index)
        dl.loc[:,'variable'] = dl.variable_o.apply(lambda x : label_dict[x])

        return [dp,dl]
        
#    
#    # Customize the exception behaviour    
    except:
        #exc.args += (cutoff)
        print('partition failed!')
        return -1
    



