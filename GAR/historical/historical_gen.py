# -*- coding: utf-8 -*-
"""
Created on Tue Jul  31 14:40:30 2018

@author: CWang2
"""

import pandas as pd          
import numpy as np
from GAR.globals import show_message
from .tskew import tskew_pdf
from .tskew import tskew_cdf
from .tskew import tskew_ppf
from .asymt import asymt_pdf
from .asymt import asymt_cdf
from .asymt import asymt_ppf
from .tskewfit import tskew_fit
from .asymtfit import asymt_fit
import matplotlib.pyplot as plt                       ## Plotting
from matplotlib.ticker import FormatStrFormatter

def historical_gen(cond_quants,fitparam,dates,realvalues,olsmeans):
    
    
    # TODO: get freq of data directly from data

    n=len(cond_quants)
    
    n_charts=10
    draws=list(range(n))
    if n>n_charts:
        draws=[int(n*i/n_charts) for i in range(n_charts)]
    print(draws)
    
    if fitparam['fittype']=='T-skew':
        tsfits=[]
        for cond_quant in cond_quants:
            tsfits.append(tskew_fit(cond_quant,fitparam))
        
        res={}
        res['location']=[]
        res['scale']=[]
        res['skew']=[]
        res['P_growth_under_0']=[]
        res['var10%']=[]
        res['var5%']=[]
        res['var50%']=[]
        res['var90%']=[]
        res['var95%']=[]
        var5=0
        var10=0
        var50=0
        var90=1
        var95=1
        chartpacks=[]
        ct=0
        for tsfit in tsfits:
            res['location'].append(tsfit['loc']/tsfit['scale'])
            res['scale'].append(tsfit['scale'])
            res['skew'].append(tsfit['skew'])
            min_v = res['location'][-1]-10
            max_v = res['location'][-1]+10
            x_list = [x for x in np.arange(min_v,max_v,0.05)]
            growth0=tskew_cdf(0, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            res['P_growth_under_0'].append(growth0)
            #ycdf = [tskew_cdf(z, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) for z in x_list]
            var5=tskew_ppf(0.05, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            var10=tskew_ppf(0.1, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            var50=tskew_ppf(0.5, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            var90=tskew_ppf(0.9, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            var95=tskew_ppf(0.95, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])            
            res['var5%'].append(var5)
            res['var10%'].append(var10)
            res['var50%'].append(var50)
            res['var90%'].append(var90)
            res['var95%'].append(var95)
            if ct in draws:
                figchart, ax = plt.subplots(1, 1, figsize=(10,5))
                yvals= [tskew_pdf(z, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) for z in x_list] 
            
                titlestr = " T-skew quantile fit for "+str(dates[ct])[:10]+" growth rate"
                lablestr = "Density "+str(dates[ct])[:10]+" "+"growth rate"
                ax.plot(x_list,yvals,'b-',label=lablestr)
                if np.isnan(realvalues[ct]):
                    modx=res['location'][-1]
                else:
                    modx=realvalues[ct]
                mody=max(yvals)
                ax.plot([modx,modx],[0,mody],'r-.')
                ax.set_title(titlestr)
                ax.legend()
                chartpacks.append(figchart)
            ct+=1
        print(res['var10%'])
        para=['var10%','location','scale','skew','P_growth_under_0']
        fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(13,80))
        
        for i,k in enumerate (para):
            axes[i].plot(dates,res[k])
            axes[i].set_title('Historical distribution of T-skew parameter: {}'.format(k), fontsize=30, y=1.05)
            axes[i].tick_params(labelsize=25)
            axes[i].set_xlabel('')
            
        return fig,res,chartpacks
    elif fitparam['fittype']=='Asymmetric T':
        asfits=[]
        
        for ind,cond_quant in enumerate(cond_quants):
            asfits.append(asymt_fit(cond_quant,fitparam,olsmeans[ind]))
        
        res={}
        res['location']=[]
        res['scale']=[]
        res['skew']=[]
        res['P_growth_under_0']=[]
        res['var10%']=[]
        res['var5%']=[]
        res['var50%']=[]
        res['var90%']=[]
        res['var95%']=[]
        var5=0
        var10=0
        var50=0
        var90=1
        var95=1
        chartpacks=[]
        ct=0
        for asfit in asfits:
            res['location'].append(asfit['loc'])
            res['scale'].append(asfit['scale'])
            res['skew'].append(asfit['skew'])
            min_v = res['location'][-1]-1.5
            max_v = res['location'][-1]+1.5
            x_list = [x for x in np.arange(min_v,max_v,0.02)]
            growth0=asymt_cdf(0, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            res['P_growth_under_0'].append(growth0)
            #ycdf = [asymt_cdf(z, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale']) for z in x_list]
            var5=asymt_ppf(0.05, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            var10=asymt_ppf(0.1, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            var50=asymt_ppf(0.5, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            var90=asymt_ppf(0.9, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            var95=asymt_ppf(0.95, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])            
            res['var5%'].append(var5)
            res['var10%'].append(var10)
            res['var50%'].append(var50)
            res['var90%'].append(var90)
            res['var95%'].append(var95)
            if ct in draws:
                figchart, ax = plt.subplots(1, 1, figsize=(10,5))
                yvals= [asymt_pdf(z, alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])for z in x_list] 
            
                titlestr = " Asymmetric T quantile fit for "+str(dates[ct])[:10]+" growth rate"
                lablestr = "Density "+str(dates[ct])[:10]+" "+"growth rate"
                ax.plot(x_list,yvals,'b-',label=lablestr)
                if np.isnan(realvalues[ct]):
                    modx=res['location'][-1]
                else:
                    modx=realvalues[ct]
                mody=max(yvals)
                ax.plot([modx,modx],[0,mody],'r-.')
                ax.set_title(titlestr)
                ax.legend()
                chartpacks.append(figchart)
            ct+=1
        print(res['var10%'])
        para=['var10%','location','scale','skew','P_growth_under_0']
        fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(20,60))
        
        for i,k in enumerate (para):
            axes[i].plot(dates,res[k])
            axes[i].set_title('Historical distribution of asymmetric T parameter: {}'.format(k), fontsize=24, y=1.05)
            axes[i].tick_params(labelsize=24)
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            axes[i].set_xlabel('')
            
        return fig,res,chartpacks