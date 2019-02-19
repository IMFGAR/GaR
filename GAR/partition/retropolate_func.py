# -*- coding: utf-8 -*-
"""
Run the quantile regressions
rlafarguette@imf.org
Time-stamp: "2018-01-11 14:39:36 rlafarguette"
"""

###############################################################################
#%% Modules
###############################################################################
## Core modules
import pandas as pd          
import numpy as np
from datetime import datetime as date      

## Self-defined modules
#import paths
#import globalvars as gv; importlib.reload(gv) ## Partitioning


# Zscore correction
def zscore(series):
    return((series - series.mean())/series.std(ddof=0))

###############################################################################
# Given two frame of signle country retroplate late frame to early frame, 
# return the retroplated frame
def retropolate(dfearly,dflate,complete_early,groups_dict):
###############################################################################

    ###########################
    ###TODO Remove country#####
    dfearly['country']=0
    dflate['country']=0
    ###########################
    
    #dload = pd.read_excel(gv.final_data_dir + '/Partitions_late.xlsx',
    #                    sheetname='Loadings') 

    ## Select the data of interest
    ## This part can be removed as it shoud be done outside of the function.
    group_vars = [x for x in groups_dict.keys()]
    all_vars = ['country', 'date'] + group_vars

    de = dfearly.loc[:,all_vars].copy()
    dl = dflate.loc[:,all_vars].copy()
## Sort the frames
    de = dfearly.sort_values(by=['country', 'date'], ascending=[1,1])
    dl = dflate.sort_values(by=['country', 'date'], ascending=[1,1])

###############################################################################
#%% For every country, compute the reverse growth rate based on early data
###############################################################################
## Compute the reverse delta (from future to now, data inverted)
    for pvar in group_vars:
        rgr_n = '{}_rgr'.format(pvar)

    ## Need to normalize: compute the zscore, per country 
        de[pvar] = de.groupby(['country'])[pvar].apply(zscore)   
        dl[pvar] = dl.groupby(['country'])[pvar].apply(zscore)
    
    ## Compute the delta, per country (pay attention to the order, future second)
        de[rgr_n] = de.groupby(['country'])[pvar].apply(lambda x: x - x.shift(-1))
        dl[rgr_n] = dl.groupby(['country'])[pvar].apply(lambda x: x - x.shift(-1))     

###############################################################################
#%% Index creation using the reverse delta
## Dulani's trick: sum for small numbers, growth rate for large number !!
###############################################################################

    # 1. Identify the missing dates from the late frame
    # dec = de.loc[de.country==pays,:]
    # dmc = dm.loc[dm.country==pays,:]
    # dlc = dl.loc[dl.country==pays,:]

    ####### From late to middle
    late_missing_dates = sorted(list(set(de.date) - set(dl.date)))
    late_start_date = min(dl.date)
    

    ## Isolate the middle frame without long time frame
    ef = de.loc[de.date < late_start_date, :]
    ef = ef.sort_values(by='date', ascending=0) # Reverse cum sum !!


    # 2. Compute the cumulative growth rate based only on the recent frame
    for pv in group_vars:
        ef['{}_cum_rgr'.format(pv)] = ef['{}_rgr'.format(pv)].cumsum().copy()

    ## 3. Using cumulative sum, create the missing frame
    mgr_frames_list = list()

    ## Retroplating for every group, only incomplete group will be updated.
    for group in group_vars:
        
        start_val = dl.loc[dl['date'] == min(dl['date']),group].values[0]
        
#        print('start_val',start_val, min(dl.date),group)
#        print(dl.iloc[0])
#        print(de.iloc[0])
        dng = pd.DataFrame(index=late_missing_dates, columns=['date'])
        dng['date'] = dng.index
        dng['country'] = de['country'].values[0]
        dng = dng.sort_values(by='date', ascending=0)

        gr_cum = '{}_cum_rgr'.format(group)
        dng_f = dng.merge(ef[['date', gr_cum]], on=['date'], how='left')
        #dng_f[group] = dng_f[gr_cum] + start_val # Increment the value

        #If group in the early frame is complete, no retroplation is needed.
        #Use the value in early group

        if group in complete_early:
            dng_f[group]=ef[group].values
#            print('complete')
#            print(dng_f[group].values)
#            print(ef[group].values)
        else:
            dng_f[group] = dng_f[gr_cum] + start_val


        dng_f.index=dng_f['date']
#        print(group)
#        print('dng_f')
#        print(dng_f)
#        print('ef')
#        print(ef)
        mgr_frames_list.append(dng_f)
        
      
    ## Merge the new groups into a early augmented frame
    dea = mgr_frames_list[0]
    for frame in mgr_frames_list[1:]:
        dea = pd.merge(dea, frame, on=['date', 'country'])
    dea.index=dea['date']

    ## Merge late, early augmented

    d_complete = pd.concat([dl[all_vars], dea[all_vars]],axis='index')
    

    d_complete=d_complete.sort_values(by=['country', 'date'])
    dfearly=dfearly.sort_values(by=['country', 'date'])
    
    ## complete group fix
    for group in complete_early:
        d_complete[group]=dfearly[group]

    return d_complete



















