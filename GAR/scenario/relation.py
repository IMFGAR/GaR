# -*- coding: utf-8 -*-
"""
Created on Tue Jul  31 14:40:30 2018

@author: CWang2
"""

import pandas as pd          
import numpy as np
from GAR.globals import show_message


# Calculate shock relations
def gen_relation(shockdict,partition_groups,df_data,df_partition):
    
    df_shockedvar = pd.DataFrame(index=df_data.index)
    df_shockedgrp = df_partition.copy()
    for group in df_shockedgrp.columns:
        df_shockedgrp[group+'_shocked']=df_shockedgrp[group]
    for var,shock in shockdict.items():
        ct=0

        if shock['shocktype']=='By +/- STD':
            if var in partition_groups.keys():
                std=np.nanstd(df_shockedgrp[var])
            else:
                df_shockedvar[var]=df_data[var]
                std=np.nanstd(df_data[var].values)
                df_shockedvar[var+'_shocked']=df_shockedvar[var]+std*shock['shockvalue']
        elif shock['shocktype']=='By +/- percentage' and (var not in partition_groups.keys()):
            df_shockedvar[var]=df_data[var]
            df_shockedvar[var+'_shocked']=df_shockedvar[var]*(1+shock['shockvalue'])
        for group,compvars in partition_groups.items():
            if var in compvars and var in partition_groups.keys():
                print(var+' is not well  defined.')
                message = var+' is in partition groups and also a group name. Please Check'
                show_message(message,halt = False) 
                
            if var in compvars:
                ct+=1
                df_var=df_data[['date',var]].dropna()
                df_part=df_partition[['date',group]].dropna()

                sdate=max(min(df_var['date'].values),(min(df_part['date'].values)))
                edate=min(max(df_var['date'].values),(max(df_part['date'].values)))

                df_var=df_var[(df_var['date']>=sdate) & (df_var['date']<=edate)]
                df_part=df_part[(df_part['date']>=sdate) & (df_part['date']<=edate)]
                
                cov=np.corrcoef(df_var[var].values,df_part[group].values)[0][1]                
                #cov=np.cov(df_data[var].values,df_partition[group])[0][1]
                print(group,cov,var,shock['shocktype'],shock['shockvalue'])
                if shock['shocktype']=='By +/- STD':
                    df_shockedgrp[group+'_shocked']= df_shockedgrp[group+'_shocked']+std*shock['shockvalue']*cov
                elif shock['shocktype']=='By +/- percentage':                
                    df_shockedgrp[group+'_shocked']= df_shockedgrp[group+'_shocked']+df_data[var]*shock['shockvalue']*cov                
            elif var==group:
                ct+=1
                print(group,var,shock['shocktype'],shock['shockvalue'])
                if shock['shocktype']=='By +/- STD':
                    df_shockedgrp[group+'_shocked']= df_shockedgrp[group+'_shocked']+std*shock['shockvalue']
                elif shock['shocktype']=='By +/- percentage':                
                    df_shockedgrp[group+'_shocked']= df_shockedgrp[group+'_shocked']+df_shockedgrp[group]*shock['shockvalue']

        if ct==0:
            print(var+' not in any group.')
            message = var+' not in any partition groups. Please Check'
            show_message(message,halt = False)   
    #print(df_shockedgrp.head())                   
    return df_shockedvar,df_shockedgrp
