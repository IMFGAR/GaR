# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:40:33 2018

@author: CWang2
"""
import pandas as pd
import numpy as np
from .tskew import tskew_pdf
from .tskew import tskew_cdf
from .tskew import tskew_ppf
from .tskew import tskew_mean
from .asymt import asymt_pdf
from .asymt import asymt_cdf
from .asymt import asymt_ppf
from .asymt import asymt_mean
from .tskewfit import tskew_fit
from .asymtfit import asymt_fit
import matplotlib.pyplot as plt  
from matplotlib.ticker import FormatStrFormatter

def gen_seg_skewt(fitdates,fitparam,skewtlist,horizonlist,medianlist,loclist,inputdatelist,freq='Quarterly'):
    
    n=len(skewtlist)
    colorlist=['red','blue','green','cyan','magenta','orange','lime','violet','crimson']
    ymax=-1
    if fitparam['fittype']=='T-skew':
        min_v = min(loclist)-8
        max_v = max(loclist)+8
        x_list = [x for x in np.arange(min_v,max_v,0.05)]
        titlestr = freq+" T-skew forecast for growth rate"
        fig, ax = plt.subplots(1, 1, figsize=(20,10))
        ax.set_title(titlestr,fontsize=24)
        ax.legend(fontsize=24)
        ax.tick_params(labelsize=24)
        plt.legend(loc=2)
        plt.ylabel('Probability Density', fontsize=24)
        plt.xlabel(freq+' GDP(compound annual growth rate)', fontsize=24)    
        
        titlestr_cdf = freq+" T-skew forecast for growth rate"
        fig2, ax2 = plt.subplots(1, 1, figsize=(20,10))
        ax2.set_title(titlestr_cdf,fontsize=24)
        ax2.set_ylim(0, 1)
        ax2.set_title(titlestr,fontsize=24)        
        ax2.legend(fontsize=24,loc=2)
        ax2.tick_params(labelsize=24)
        plt.legend(loc=2)
        plt.xlim(x_list[0],x_list[-1])
        plt.ylabel('Cumulative probability', fontsize=24)
        plt.xlabel(freq+' GDP(compound annual growth rate)', fontsize=24)
        df_header=['Tskew_PDF_x']
    elif fitparam['fittype']=='Asymmetric T':
        min_v = min(loclist)-1.5
        max_v = max(loclist)+1.5
        x_list = [x for x in np.arange(min_v,max_v,0.01)]
        titlestr = freq+" Asymmetric T forecast for growth rate"
        fig, ax = plt.subplots(1, 1, figsize=(20,10))
        ax.set_title(titlestr,fontsize=24)
        ax.legend(fontsize=24)
        ax.tick_params(labelsize=24)
        plt.legend(loc=2)
        plt.ylabel('Probability Density', fontsize=24)
        plt.xlabel(freq+' GDP(compound annual growth rate)', fontsize=24)    
        
        titlestr_cdf = freq+" Asymmetric T forecast for growth rate"
        fig2, ax2 = plt.subplots(1, 1, figsize=(20,10))
        ax2.set_title(titlestr_cdf,fontsize=24)
        ax2.set_ylim(0, 1)
        ax2.set_title(titlestr,fontsize=24)        
        ax2.legend(fontsize=24,loc=2)
        ax2.tick_params(labelsize=24)
        plt.legend(loc=2)
        plt.xlim(x_list[0],x_list[-1])
        plt.ylabel('Cumulative probability', fontsize=24)
        plt.xlabel(freq+' GDP(compound annual growth rate)', fontsize=24)
        df_header=['AsymT_PDF_x']
    
    
    df_tmp=[x_list]
    horizons=['Forward horizon']
    horizons.extend(horizonlist)
    inputdate=['Input data cut off']
    inputdate.extend(inputdatelist)
    cmode=['Conditional mode']
    cmedian=['Conditional median']
    cmean=['Conditional mean']
    gar5=['GaR5%']
    gar10=['GaR10%']
    gzero=['Growth below 0 probablity']
    xq5s=[]
    for indhz in range(n):
        
        if fitparam['fittype']=='T-skew':
            tsfit=skewtlist[indhz]
            yvals= [tskew_pdf(z, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) for z in x_list]    
            ymax=max(ymax,max(yvals))
            ycdf = [tskew_cdf(z, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) for z in x_list]
            yzero=tskew_cdf(0, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            ax.plot(x_list,yvals,'-',color=colorlist[indhz],label=fitdates[indhz].strftime('%m/%d/%Y')+" forward "+str(horizonlist[indhz]))
            ax2.plot(x_list,ycdf,'-',color=colorlist[indhz],label=fitdates[indhz].strftime('%m/%d/%Y')+" forward "+str(horizonlist[indhz]))
            df_header.append('Tskew_PDF_y_PROJ'+str(indhz+1))
            df_header.append('Tskew_CDF_y_PROJ'+str(indhz+1))
            df_tmp.append(yvals)
            df_tmp.append(ycdf)
        
            xq5=tskew_ppf(0.05, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) 
            xq10=tskew_ppf(0.1, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) 
            yq5= tskew_pdf(xq5, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) 
            yq10= tskew_pdf(xq10, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])     
            ycq5= tskew_cdf(xq5, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew']) 
            ycq10= tskew_cdf(xq10, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            meanx=tskew_mean(df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            modx=tsfit['loc']
            medx=tskew_ppf(0.5, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
            xq5s.append(xq10)
    
        elif fitparam['fittype']=='Asymmetric T':
            asymtfit=skewtlist[indhz]
            yvals= [asymt_pdf(z, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) for z in x_list]    
            ymax=max(ymax,1.2*max(yvals))
            ycdf = [asymt_cdf(z, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) for z in x_list]
            
            ax.plot(x_list,yvals,'-',color=colorlist[indhz],label=fitdates[indhz].strftime('%m/%d/%Y')+" forward "+str(horizonlist[indhz]))
            ax2.plot(x_list,ycdf,'-',color=colorlist[indhz],label=fitdates[indhz].strftime('%m/%d/%Y')+" forward "+str(horizonlist[indhz]))
            df_header.append('AsymT_PDF_y_PROJ'+str(indhz+1))
            df_header.append('AsymT_CDF_y_PROJ'+str(indhz+1))
            df_tmp.append(yvals)
            df_tmp.append(ycdf)
            xq5 = asymt_ppf(0.05, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) 
            xq10 = asymt_ppf(0.1, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale'])     
            yq5= asymt_pdf(xq5, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) 
            yq10= asymt_pdf(xq10,alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) 
            ycq5= asymt_cdf(xq5, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) 
            ycq10= asymt_cdf(xq10,alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) 
            yzero=asymt_cdf(0,alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale']) 
            meanx=asymt_mean(alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale'])
            modx=loclist[indhz]
            medx=asymt_ppf(0.5, alpha=asymtfit['skew'], nu1=asymtfit['kleft'], nu2=asymtfit['kright'], mu=asymtfit['loc'], sigma=asymtfit['scale'])
            xq5s.append(xq5)
            
        cmode.append(float("{:.4f}".format(modx)))
        cmedian.append(float("{:.4f}".format(medx)))
        cmean.append(float("{:.4f}".format(meanx)))
        gar5.append(float("{:.4f}".format(xq5)))
        gar10.append(float("{:.4f}".format(xq10)))
        gzero.append(float("{:.4f}".format(yzero)))
            
            
    dfpdf=pd.DataFrame(df_tmp)
    dfpdf = dfpdf.transpose()
    dfpdf.columns = df_header
    if fitparam['fittype']=='T-skew':
        ax.set_ylim(0, 1.2*ymax)
    elif fitparam['fittype']=='Asymmetric T':
        ax.set_ylim(0,1.05*ymax)
        c=(min(loclist)+max(loclist))/2
        l=max(xq5s)
        ax.set_xlim(l,1.8*c-0.8*l)
        ax2.set_xlim(l,1.8*c-0.8*l)
    ax.legend(fontsize=24,loc=2)
    ax2.legend(fontsize=24,loc=2)
    res=[horizons,inputdate,cmode,cmedian,cmean,gar5,gar10,gzero]
    return res,fig,fig2,dfpdf
   