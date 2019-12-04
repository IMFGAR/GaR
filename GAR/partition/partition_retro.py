# -*- coding: utf-8 -*-
"""
Partition the data for Peru using cutoffs
rlafarguette@imf.org
Time-stamp: "2018-01-08 19:38:27 rlafarguette"
Edited by cwang2@imf.org for retroplating partition.
"""

###############################################################################
#%% Modules
###############################################################################
## Core modules
import os, sys, importlib
import pandas as pd          
import numpy as np

## Self-defined modules
#import paths
#import globalvars as gv; importlib.reload(gv) ## Partitioning

from .gen_cutoff import gen_cutoff
from .retropolate_func import retropolate
from .partition_cutoff import  p_cutoff
from .cpdgrowth import cum_gr
from .cpdgrowth import yoy_gr
## Functions loading
#from partition import Partition
from datetime import datetime as date                 ## Dates
#from dateutil.relativedelta import relativedelta      ## Dates computations

## Suppress warnings
import warnings
warnings.filterwarnings("ignore")

# Zscore correction
def zscore(series):
    return((series - series.mean())/series.std(ddof=0))

###############################################################################
#Function to generate retropolated partition in a time period
###############################################################################
def partition_retro(**kwargs):

    #dall="NA",groups_dict={},tdep="NA",bench='NA',rgdp='NA', horizon=4, method='LDA',sdate=date(1,1,1),edate=date(9999,12,30), benchcutoff=0.30, saveim=False):
    if 'dall' in kwargs:
        dall=kwargs['dall']
    else:
        dall=        print('Error! No data imported')
    
    if 'groups_dict' in kwargs:
        groups_dict=kwargs['groups_dict']
    else:
        groups_dict={}
        print('Warning :group dict not specified')
        
    if 'tdep' in kwargs:
        tdep=kwargs['tdep']
    else:
        tdep='NA'
           
    if 'rgdp' in kwargs:
        rgdp=kwargs['rgdp']
    else:
        rgdp='NA'
        print('Warning :real gdp not specified')
    
    if 'method_growth' in kwargs:
        method_growth=kwargs['method_growth']
    else:
        method_growth='cpd' #'yoy'
        
    if 'horizon' in kwargs:
        horizon=kwargs['horizon']
    else:
        horizon=4
        print('Warning :horizon not specified')
    
    if 'method' in kwargs:
        method=kwargs['method']
    else:
        method='LDA'
        print('Warning :method not specified, using LDA')
        
    if 'sdate' in kwargs:
        sdate=kwargs['sdate']
    else:
        sdate=date(1,1,1)
        print('Warning : start date specified')
    
    if 'edate' in kwargs:
        edate=kwargs['edate']
    else:
        edate=date(9999,12,30)
        print('Warning : end date specified')
        
    if 'benchcutoff' in kwargs:
        benchcutoff=kwargs['benchcutoff']
    else:
        benchcutoff=0.3        
        print('Warning : bench mark cutoff not specified')
        
    if 'PLStarget' in kwargs:
        PLStarget=kwargs['PLStarget']
    else:
        PLStarget=None
        
    if 'saveim' in kwargs:
        saveim=kwargs['saveim']
    else:
        saveim=False        
    

        
## Some data treatment for peru (I have some missing data at the end...)
    log_frame=pd.DataFrame(columns=['Time','Action'])
    dall = dall.fillna(method='ffill').copy()
    dall = dall[(dall['date']>=sdate) & (dall['date']<=edate)]
    
## Calculating the growth
    dall['dummygrp']=1    
    if method_growth=='cpd':
        dall.loc[:,tdep] = dall.groupby(['dummygrp'])[rgdp].apply(cum_gr,horizon=horizon)
    elif method_growth=='yoy':
        dall.loc[:,tdep] = dall.groupby(['dummygrp'])[rgdp].apply(yoy_gr,horizon=horizon)
        if horizon<4:
            dall = dall.iloc[4-horizon-1:]
            sdate=dall.index.values[0]
    elif method_growth=='level':
        dall.loc[:,tdep] = dall.groupby(['dummygrp'])[rgdp].shift(-horizon)
    else:
        # Assume data provided in the data sheet
        pass
    
    

## Using dependent variable as benchmark, although it is an extra copy, the code for benchmark is
## written for accepting any bench mark variable, and keeping this will give the flexibility for
## futrue using other variables.
    bench='benchvar'
    dall[bench]=dall[tdep]

## Generating all cutoffs in the period, sorted from latest to earliest
    [cutoffs,complete_group]=gen_cutoff(dall=dall,groups_dict=groups_dict,startdate=sdate,enddate=edate)

    if (cutoffs==-1):
        tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
        action="In the given time period some groups are complete empty. No feasible partition can be made."
        log = pd.Series({'Time': tn, 'Action': action})
        log_frame=log_frame.append(log,ignore_index=True)
        return dall.head(),dall.head(),log_frame,-1

    if len(cutoffs)==0:
        tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
        action="No data in the cutoff period"
        log = pd.Series({'Time': tn, 'Action': action})
        log_frame=log_frame.append(log,ignore_index=True)
        return dall.head(),dall.head(),log_frame,-1

    
## Generating the parition for the latest cutoff            
    [dp1,dl]=p_cutoff(dall,groups_dict,cutoffs[0],bench,method, benchcutoff, PLStarget, saveim=saveim)

#    raw_name = r'\Retropolated_partitions_'+country+'.xlsx'
#    xlname = gv.final_data_dir + raw_name
    
## If only one cutoff, no retroplating needed, Just save the result.
#    if len(cutoffs)==1:       
#        ## Save the results
##        writer = pd.ExcelWriter(xlname, engine='xlsxwriter')
##        dp1.to_excel(writer, 'Partition data', index=False)
##        dl.to_excel(writer, 'Loadings', index=False) # From the latest frame
##        writer.save()
#        tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
#        action="Data are complete in the time periord. No need for retroplation." 
#        log = pd.Series({'Time': tn, 'Action': action})
#        log_frame=log_frame.append(log,ignore_index=True)
#        return dp1,dl,log_frame,1

#        print('Partition results saved')
        
## Repeatedly generate partition of new cutoff, retropolate the previous partition to it.

    dpo=dp1

    
    for i in range(1,len(cutoffs)):

        [dpn,dln]=p_cutoff(dall,groups_dict,cutoffs[i],bench,method, benchcutoff, PLStarget, saveim=saveim)    

        dpr=retropolate(dfearly=dpn,dflate=dpo,complete_early=complete_group[i],groups_dict=groups_dict)


        dpo=dpr.copy()
        tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
        retrovar=" "
        for e in groups_dict:
            if e not in complete_group[i]:
                retrovar+=e+", "
        action="Retroplating for "  + retrovar
        log = pd.Series({'Time': tn, 'Action': action})
        log_frame=log_frame.append(log,ignore_index=True)

    dl['cutoff']=sdate
    if method=='PLS':
        dl=dl[['variable','cutoff','loadings','group','vip']]
    else:
        dl=dl[['variable','cutoff','loadings','group','variance_ratio']]
        
## Compute the zscore for the final frame to makes them consistent
    group_vars = [x for x in groups_dict.keys()]
    for group in group_vars:
        #print(group)
        #print(dpo[group].iloc[0:10])
        v=dpo[group].values
        print('before',group,np.var(v))
        dpo[group] = zscore(dpo[group])
        #print(dpo[group].iloc[0:10])
        v=dpo[group].values
        print('after',group,np.var(v))
        
    dpo.index.name=None
    dall.index.name=None
    dretro_final = dpo .merge(dall[['date',tdep]], on=['date'], how='left')
    dretro_final.index=dretro_final['date']
    dretro_final.index.name=None

        
    tn=date.now().strftime('%Y-%m-%d %H:%M:%S')
    action="Retroplating successfully finished." 
    log = pd.Series({'Time': tn, 'Action': action})
    log_frame=log_frame.append(log,ignore_index=True)
#    writer = pd.ExcelWriter(xlname, engine='xlsxwriter')
#    dretro_final.to_excel(writer, 'Partition data', index=False)
#    dl.to_excel(writer, 'Loadings', index=False) # From the latest frame   
#    print("Retroplating partition for "+country+" saved!")
#    writer.save()
#    existcode=1
    #print(dretro_final.head())
    #print(dretro_final.columns)
    #print(dretro_final)
    try:
        dretro_final.drop(['country'],axis=1,inplace=True)
    except:
        pass
    return dretro_final,dl,log_frame,1
#
#if __name__ == "__main__":
#    partition_retro("PER","NA",groups_dict,date(1995,1,1),date(2017,8,30),saveim=False)
