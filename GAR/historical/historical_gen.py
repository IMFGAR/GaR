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
from math import log

def historical_gen(cond_quants,fitparam,dates,realvalues,olsmeans):
    
    
    # TODO: get freq of data directly from data

    n=len(cond_quants)
    
    n_charts=10
    draws=list(range(n))
    if n>n_charts:
        draws=[int(n*i/n_charts) for i in range(n_charts)]

    
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
        pits=[]
        logscore={}
        logscore['uncensored']=[]
        logscore['10tail']=[]
        logscore['90tail']=[]
        
        for ct,tsfit in enumerate(tsfits):
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
            
            if not np.isnan(realvalues[ct]):
                realcdf=tskew_cdf(realvalues[ct], df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
                realpdf=tskew_pdf(realvalues[ct], df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            else:
                realcdf=tskew_cdf(tsfit['loc']/tsfit['scale'], df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
                realpdf=tskew_pdf(tsfit['loc']/tsfit['scale'], df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            pits.append(realcdf)
            logscore['uncensored'].append(log(realpdf))
            if realcdf<=0.1:
                logscore['10tail'].append(log(realpdf))
            else:
                logscore['10tail'].append(log(0.9))
            if realcdf>=0.9:
                logscore['90tail'].append(log(realpdf))
            else:
                logscore['90tail'].append(log(0.9))

                
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
        
        pits_cdf=[]
        npits=len(pits)
        for r in np.arange(0,1,0.01):
            pits_cdf.append(len([x for x in pits if x<=r])/npits) # Calculate how many realized values are below any given probability
        
        figpit, axpit= plt.subplots(1, 1, figsize=(8,8))
        axpit.plot(list(np.arange(0,1,0.01)),pits_cdf,'r-',label='Realized')
        axpit.plot(list(np.arange(0,1,0.01)),list(np.arange(0,1,0.01)),'b-',label='U~(0,1)')
        axpit.plot(list(np.arange(0,1,0.01)),[e+1.34*npits**(-0.5) for e in np.arange(0,1,0.01)],'b-',label='5 percent critical values',linestyle='dashed')
        axpit.plot(list(np.arange(0,1,0.01)),[e-1.34*npits**(-0.5) for e in np.arange(0,1,0.01)],'b-',linestyle='dashed')
        axpit.set_title('Probablity inversion test', fontsize=16, y=1.01)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.legend(loc=4,fontsize=12)

        figls, axls= plt.subplots(1, 1, figsize=(12,8))
        axls.plot(dates, logscore['uncensored'], 'k-', label='Uncensored log score')
        axls.plot(dates, logscore['10tail'], 'r--',label='Censored log score for left  10% tail')
        axls.plot(dates, logscore['90tail'], 'g-.', label='Censored log score for right 10% tail')
        plt.legend(loc=3,fontsize=12)
        plt.title('Uncensored and Censored Logscores',fontsize=16,y=1.01)
        
        
        para=['var10%','location','scale','skew','P_growth_under_0']
        fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(15,80))
        
        for i,k in enumerate (para):
            axes[i].plot(dates,res[k])
            axes[i].set_title('Historical distribution of T-skew parameter: {}'.format(k), fontsize=24, y=1.03)
            axes[i].tick_params(labelsize=16)
            axes[i].set_xlabel('')
        figs={}   
        figs['res']=fig
        figs['pit']=figpit
        figs['ls']=figls
        return figs,res,chartpacks
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
        pits=[]
        logscore={}
        logscore['uncensored']=[]
        logscore['10tail']=[]
        logscore['90tail']=[]
        for ct,asfit in enumerate(asfits):
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
            if not np.isnan(realvalues[ct]):
                realcdf=asymt_cdf(realvalues[ct], alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
                realpdf=asymt_pdf(realvalues[ct],alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            else:
                realcdf=asymt_cdf(asfit['loc'], alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
                realpdf=asymt_pdf(asfit['loc'], alpha=asfit['skew'], nu1=asfit['kleft'], nu2=asfit['kright'], mu=asfit['loc'], sigma=asfit['scale'])
            pits.append(realcdf)
            logscore['uncensored'].append(log(realpdf))
            if realcdf<=0.1:
                logscore['10tail'].append(log(realpdf))
            else:
                logscore['10tail'].append(log(0.9))
            if realcdf>=0.9:
                logscore['90tail'].append(log(realpdf))
            else:
                logscore['90tail'].append(log(0.9))
                
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
                
        pits_cdf=[]        
        npits=len(pits)
        for r in np.arange(0,1,0.01):
            pits_cdf.append(len([x for x in pits if x<=r])/npits) # Calculate how many realized values are below any given probability
        
        figpit, axpit= plt.subplots(1, 1, figsize=(8,8))
        axpit.plot(list(np.arange(0,1,0.01)),pits_cdf,'r-',label='Realized')
        axpit.plot(list(np.arange(0,1,0.01)),list(np.arange(0,1,0.01)),'b-',label='U~(0,1)')
        axpit.plot(list(np.arange(0,1,0.01)),[e+1.34*npits**(-0.5) for e in np.arange(0,1,0.01)],'b-',label='5 percent critical values',linestyle='dashed')
        axpit.plot(list(np.arange(0,1,0.01)),[e-1.34*npits**(-0.5) for e in np.arange(0,1,0.01)],'b-',linestyle='dashed')
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.legend(loc=4)
        plt.title('Probablity inversion test.')
        
        figls, axls= plt.subplots(1, 1, figsize=(12,8))
        axls.plot(dates, logscore['uncensored'], 'k-', label='Uncensored log score')
        axls.plot(dates, logscore['10tail'], 'r--',label='Censored log score for left  10% tail')
        axls.plot(dates, logscore['90tail'], 'g-.', label='Censored log score for right 10% tail')
        plt.legend(loc=3,fontsize=12)
        plt.title('Uncensored and Censored Logscores',fontsize=16,y=1.01)
        
        para=['var10%','location','scale','skew','P_growth_under_0']
        fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(15,80))
        
        for i,k in enumerate (para):
            axes[i].plot(dates,res[k])
            axes[i].set_title('Historical distribution of asymmetric T parameter: {}'.format(k), fontsize=24, y=1.03)
            axes[i].tick_params(labelsize=16)
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            axes[i].set_xlabel('')
        figs={}   
        figs['res']=fig
        figs['pit']=figpit   
        figs['ls']=figls
        return figs,res,chartpacks