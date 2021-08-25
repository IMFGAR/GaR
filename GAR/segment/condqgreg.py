# -*- coding: utf-8 -*-
"""
Created on Thu May  3 11:01:30 2018

@author: cwang2
"""

from datetime import datetime as date

import pandas as pd
import numpy as np
from sklearn.preprocessing import scale               ## Scaling (zscore)

from .quantilereg import QuantileReg

def condquant(dall,depvar,regressors_avl,horizon,ql):
#if 1==1:    
    ql.sort()
    c_id_dict = {'horizon' : horizon}
#    orgdep=depvar    
#    dpv=list(dall[orgdep].values)
#    dpv=dpv[horizon:]
#    ev=dpv[-1]
#    for i in range(horizon):
#        dpv.append(ev)    
#    depvar=depvar+'_cpd_'+str(horizon)
#    dall[depvar]=dpv    
#    print(dall[depvar].values[:10])
    
    dall=dall.dropna(subset=regressors_avl)
    print(regressors_avl)
    #dall=dall.dropna()

    qrs = QuantileReg(depvar, indvars=regressors_avl,
                      quantile_list=ql,
                      data=dall,
                      scaling=True, alpha=0.1)

    dc = qrs.coeff
    dc = add_id(dc,c_id_dict)
    dc.insert(0, 'variable', dc.index)
    
        ## Without scaling: get the conditional quantiles 
    qru = QuantileReg(depvar, indvars=regressors_avl,
                      quantile_list=ql,
                      data=dall,
                      scaling=False, alpha=0.1)

        ## Run the predictions on the full frame (estimates can differ)

    dcq = qru.cond_quantiles(predictors = dall).copy()
#    dcq['date'] = dall['date']
#    dcq.index=dcq['date']
    dcq = add_id(dcq, c_id_dict)

    print('Quantile reg done')
        ## Store the coefficients
    dci = qru.coeff; dci = add_id(dci, c_id_dict)
    dci.insert(0, 'variable', dci.index)
    dc.rename(columns={'coeff':'coeff_scale'},inplace=True)
    dc['coeff_noscale']=dci['coeff']
    dc=dc[['variable','horizon','quantile','coeff_scale','coeff_noscale','pval','lower','upper','R2_in_sample']]

    exitcode=1
    return [dc,exitcode]

def add_id(df, id_dict):
    """ Add identifiers variables to a pandas frame """
    variables_id = sorted(list(id_dict.keys()))
    for v, var in enumerate(variables_id):
        df.insert(v, var, id_dict[var])
    return(df)

