# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 15:32:04 2018
Annualized compound GDP Growth calculation
@author: cwang2
"""

def cum_gr(series, horizon ,yearfreq=4): 
    ## Compute the compound annualized quarterly growth rate over a certain horizon
    cagr = ((series.shift(-horizon)/series)**(1/horizon))-1
    ## Need to annualize it now
    annual_cagr = ((1+cagr)**yearfreq) -1
    return(100*annual_cagr)


def yoy_gr(series, horizon, yearfreq=4): 
    ## We assume that the growth rate is quarterly. In the future, rather than having +4, should use an index period
    yoy_gr = (series.shift(-horizon)/series.shift(-horizon+yearfreq))-1
    return(100*yoy_gr)
