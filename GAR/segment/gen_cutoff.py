# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 12:47:52 2018

@author: cwang2
"""

import os, sys, importlib
import pandas as pd          
import numpy as np

## Self-defined modules
#import paths
#import globalvars as gv; importlib.reload(gv) ## Partitioning
from GAR.partition import partition

## Functions loading
from datetime import datetime as date                 ## Dates
from dateutil.relativedelta import relativedelta      ## Dates computations

## Suppress warnings
import warnings
warnings.filterwarnings("ignore")


## Function to generate partition cutoff points and completed groups
## at the coressponding cutoff ponit. Completed group will not be retropolated.
def gen_cutoff (dall="default",groups_dict={}, startdate=date(year=1,month=1,day=1), enddate=date(year=9999,month=12,day=31),):
    


    if len(dall)==0:
        print("No data found")
        return -1
    
    dall = dall.fillna(method='ffill').copy()
    dall = dall[dall.date>startdate]
    dall = dall[dall.date<enddate]
    partition_dict = groups_dict # Variables per group (price, leverage, etc.)
    
    set_date=set([])
    for key, values in partition_dict.items():
        for v in values:
            t=dall[v].first_valid_index() 
            if t not in set_date:
                set_date.add(t)
    dates=list(set_date)
    dates.sort(reverse=True)
    complete_groups=[]

    for d in dates:                
        tmp_c_key=[]
        for key, values in partition_dict.items():
            
            complete_key=True
            empty_key=True
            for v in values:
                if dall[v].first_valid_index()>d:
                    complete_key=False
                else:
                    empty_key=False
            if empty_key:
            #########################################################
            # There exists a complete emptry group, the cutoff date is not applicable
                print(d, key)
                for v in values:
                    print(v,dall[v].first_valid_index())
                print("In the given time period some groups are complete empty. No feasible partition can be made")
                return -1,-1
            else:
                if complete_key:
                    tmp_c_key.append(key)
        complete_groups.append(tmp_c_key)
            
            
    return dates,complete_groups


