# -*- coding: utf-8 -*-
"""
Created on Tue Jul  31 14:40:30 2018

@author: CWang2
"""

import pandas as pd          
import numpy as np
from .tskew import tskew_pdf
from .tskew import tskew_cdf
from .tskew import tskew_ppf
from .asymt import asymt_pdf
from .asymt import asymt_cdf
from .asymt import asymt_ppf
from .tskewfit import tskew_fit
from .asymtfit import asymt_fit
import matplotlib.pyplot as plt                       ## Plotting

def scenario_compare(cond_quant_raw,cond_quant_shocked,fitparam,fitparam_shocked,horizon,fitdate,ols_raw,ols_shocked):
    
    
    # TODO: get freq of data directly from data
    freq    = 'Quarterly'
    if fitparam['fittype']=='T-skew':
        

     

        print ('cond_quant_raw',cond_quant_raw)
        print ('cond_quant_shocked',cond_quant_shocked)
        tsfit_raw=tskew_fit(cond_quant_raw,fitparam)
        tsfit_shocked=tskew_fit(cond_quant_shocked,fitparam_shocked)
        loc_raw=tsfit_raw['loc']
        loc_shocked=tsfit_shocked['loc'] 
    
        print(tsfit_raw)
        print(tsfit_shocked)
        tsfit=tsfit_raw
        v_q5=tskew_ppf(0.05, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
        v_q40=tskew_ppf(0.4, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
        v_q60=tskew_ppf(0.6, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
        v_q95=tskew_ppf(0.95, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])

        min_v = v_q5-abs(v_q5-v_q40)
        max_v = v_q95+abs(v_q95-v_q60)
        
        
        tsfit=tsfit_shocked
        v_q5=tskew_ppf(0.05, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
        v_q40=tskew_ppf(0.4, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
        v_q60=tskew_ppf(0.6, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
        v_q95=tskew_ppf(0.95, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])

        min_v = min(min_v,v_q5-abs(v_q5-v_q40))
        max_v = max(max_v,v_q95+abs(v_q95-v_q60))
        
        x_list = [x for x in np.arange(min_v,max_v,0.05)]
        yvals_raw= [tskew_pdf(z, df=tsfit_raw['df'], loc=tsfit_raw['loc'], scale=tsfit_raw['scale'], skew=tsfit_raw['skew']) for z in x_list]
        yvals_shocked= [tskew_pdf(z, df=tsfit_shocked['df'], loc=tsfit_shocked['loc'], scale=tsfit_shocked['scale'], skew=tsfit_shocked['skew']) for z in x_list]
        ycdf_raw = [tskew_cdf(z, df=tsfit_raw['df'], loc=tsfit_raw['loc'], scale=tsfit_raw['scale'], skew=tsfit_raw['skew']) for z in x_list]
        ycdf_shocked = [tskew_cdf(z, df=tsfit_shocked['df'], loc=tsfit_shocked['loc'], scale=tsfit_shocked['scale'], skew=tsfit_shocked['skew']) for z in x_list]
        yzero=tskew_cdf(0, df=tsfit_raw['df'], loc=tsfit_raw['loc'], scale=tsfit_raw['scale'], skew=tsfit_raw['skew'])
        yzero_shocked=tskew_cdf(0, df=tsfit_shocked['df'], loc=tsfit_shocked['loc'], scale=tsfit_shocked['scale'], skew=tsfit_shocked['skew'])
        tmp_dic={'Tskew_PDF_x':x_list,'Tskew_PDF_y_before':yvals_raw,'Tskew_CDF_y_before':ycdf_raw,'Tskew_PDF_y_after':yvals_shocked,'Tskew_CDF_y_after':ycdf_shocked}
        dfpdf=pd.DataFrame(tmp_dic)[['Tskew_PDF_x','Tskew_PDF_y_before','Tskew_CDF_y_before','Tskew_PDF_y_after','Tskew_CDF_y_after']]
    #tmp_dic={'Tskew_PDF_x':x_list,'Tskew_PDF_y':yvals,'Tskew_CDF':ycdf}
    #dfpdf=pd.DataFrame(tmp_dic)
        q5loc_raw=min_v
        q10loc_raw=min_v
        q5loc_shocked=min_v
        q10loc_shocked=min_v
        for i,y in enumerate(ycdf_raw):
            if y>0.05:
                q5loc_raw=i
                break
        for i,y in enumerate(ycdf_raw):
            if y>0.1:
                q10loc_raw=i
                break
        
        for i,y in enumerate(ycdf_shocked):
            if y>0.05:
                q5loc_shocked=i
                break
        for i,y in enumerate(ycdf_shocked):
            if y>0.1:
                q10loc_shocked=i
                break     
        yq5_raw= tskew_pdf(x_list[q5loc_raw-1], df=tsfit_raw['df'], loc=tsfit_raw['loc'], scale=tsfit_raw['scale'], skew=tsfit_raw['skew']) 
        yq10_raw= tskew_pdf(x_list[q10loc_raw-1], df=tsfit_raw['df'], loc=tsfit_raw['loc'], scale=tsfit_raw['scale'], skew=tsfit_raw['skew']) 
    
        yq5_shocked= tskew_pdf(x_list[q5loc_shocked-1], df=tsfit_shocked['df'], loc=tsfit_shocked['loc'], scale=tsfit_shocked['scale'], skew=tsfit_shocked['skew']) 
        yq10_shocked= tskew_pdf(x_list[q10loc_shocked-1], df=tsfit_shocked['df'], loc=tsfit_shocked['loc'], scale=tsfit_shocked['scale'], skew=tsfit_shocked['skew']) 
    
#    for i,x in enumerate(x_list):
#        if x>5:
#            v5loc=i
#            break
        
        titlestr = "Scenario test for "+fitdate.strftime('%m/%d/%Y')+" "+"growth rate"+" forward "+str(horizon)
        lablestr_raw = "Density before shock"
        lablestr_shocked = "Density after shock"
        fig, ax = plt.subplots(1, 1, figsize=(20,10))

#    if fitparam['mode']['constraint']=='Free':
#        loc=tsfit['loc']/tsfit['scale']
#    modx=loc
        ax.set_title(titlestr,fontsize=24)

        if yvals_raw[q5loc_raw]>yvals_shocked[q5loc_raw]:
            
            ax.fill_between(x_list[:q5loc_raw], 0, yvals_raw[:q5loc_raw],  facecolor='c', interpolate=True)
            ax.fill_between(x_list[:q5loc_shocked], 0, yvals_shocked[:q5loc_shocked],  facecolor='r', interpolate=True)
        
        else:
            ax.fill_between(x_list[:q5loc_shocked], 0, yvals_shocked[:q5loc_shocked],  facecolor='r', interpolate=True)
            ax.fill_between(x_list[:q5loc_raw], 0, yvals_raw[:q5loc_raw],  facecolor='c', interpolate=True)
            
        ax.plot(x_list,yvals_raw,'c-',label=lablestr_raw)
        ax.plot(x_list,yvals_shocked,'r-',label=lablestr_shocked)

    
#    if fitparam['mode']['constraint']=='Free':
#        loc=tsfit['loc']/tsfit['scale']
#    modx=loc
#    mody=tskew_pdf(loc, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
#    ax.plot([modx,modx],[0,mody],'r-.')
#    ax.annotate('Mode', xy=(modx, mody),xycoords='data',
#                xytext=(modx-1.5, mody+0.1), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                connectionstyle="arc3"),)
#    ax.plot([x_list[q5loc_raw-1],x_list[q5loc_raw-1]],[0,yq5_raw],'k--')
#    ax.plot([x_list[q5loc_shocked-1],x_list[q5loc_shocked-1]],[0,yq5_shocked],'k--')
#    ax.annotate('GaR 5%, Before', xy=(x_list[max(0,q5loc_raw-1)], yq5_raw), xycoords='data',
#                xytext=(x_list[max(0,q5loc_raw-30)]-2, yq5_raw+0.1), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                                connectionstyle="angle3,angleA=90,angleB=0",facecolor="cyan"),)
#    
#    ax.annotate('GaR 5%, After', xy=(x_list[max(0,q5loc_shocked-1)], yq5_shocked), xycoords='data',
#                xytext=(x_list[max(0,q5loc_shocked-30)]-2, yq5_raw+0.13), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                                connectionstyle="angle3,angleA=90,angleB=0",facecolor="red"),)
#    ax.plot([x_list[q10loc-1],x_list[q10loc-1]],[0,yq10],'k--')
#    ax.annotate('GaR 10%', xy=(x_list[q10loc-1], yq10), xycoords='data',
#                xytext=(x_list[q10loc-1]-1.5, yq10+0.1), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                connectionstyle="angle3,angleA=0,angleB=90"),)
        ax.legend(fontsize=24)
        ax.tick_params(labelsize=24)
        plt.ylim(0, max(max(yvals_raw),max(yvals_shocked))*1.2)
        plt.ylabel('Probability Density', fontsize=24)
        plt.xlabel(freq+' GDP(compound annual growth rate)', fontsize=24)    
    
    
        res=[]
        res.append([' ','Before shock','After shock'])
        res.append(['Date of input',fitdate,fitdate])
        res.append(['Horizon forward',horizon,horizon])
        res.append(['Conditional mode',loc_raw,loc_shocked])
        res.append(['GaR5%',x_list[q5loc_raw-1],x_list[q5loc_shocked-1]])
        res.append(['GaR10%',x_list[q10loc_raw-1],x_list[q10loc_shocked-1]])
        res.append(['Growth below 0 probablity',float("{:.4f}".format(yzero)), float("{:.4f}".format(yzero_shocked))])
        res.append(['Skewness',tsfit_raw['skew'],tsfit_shocked['skew']])
        res.append(['Scale',tsfit_raw['scale'],tsfit_shocked['scale']])   
        return fig,res,dfpdf
    
    
    elif fitparam['fittype']=='Asymmetric T':
        if fitparam['mode']['constraint']=='Fixed':
            loc_raw=fitparam['mode']['value']
        else:
            loc_raw=cond_quant_raw [0.5]

    
        if fitparam_shocked['mode']['constraint']=='Fixed':
            loc_shocked=fitparam['mode']['value']
        else:
            loc_shocked=cond_quant_shocked [0.5]
     

        print ('cond_quant_raw',cond_quant_raw)
        print ('cond_quant_shocked',cond_quant_shocked)
        asfit_raw=asymt_fit(cond_quant_raw,fitparam,ols_raw)
        asfit_shocked=asymt_fit(cond_quant_shocked,fitparam_shocked,ols_shocked)
        loc_raw=asfit_raw['loc']
        loc_shocked=asfit_shocked['loc']

        min_v = min(loc_raw,loc_shocked)-1.5
        max_v = max(loc_raw,loc_shocked)+1.5
        while asymt_cdf(min_v+1, alpha=asfit_raw['skew'], nu1=asfit_raw['kleft'], nu2=asfit_raw['kright'], mu=asfit_raw['loc'], sigma=asfit_raw['scale'])>0.05 or asymt_cdf(min_v+1, alpha=asfit_shocked['skew'], nu1=asfit_shocked['kleft'], nu2=asfit_shocked['kright'], mu=asfit_shocked['loc'], sigma=asfit_shocked['scale'])>0.05:
            min_v-=0.5

        
        x_list = [x for x in np.arange(min_v,max_v,0.02)]
        yvals_raw= [asymt_pdf(z, alpha=asfit_raw['skew'], nu1=asfit_raw['kleft'], nu2=asfit_raw['kright'], mu=asfit_raw['loc'], sigma=asfit_raw['scale']) for z in x_list]        
        yvals_shocked= [asymt_pdf(z, alpha=asfit_shocked['skew'], nu1=asfit_shocked['kleft'], nu2=asfit_shocked['kright'], mu=asfit_shocked['loc'], sigma=asfit_shocked['scale']) for z in x_list]
        ycdf_raw= [asymt_cdf(z, alpha=asfit_raw['skew'], nu1=asfit_raw['kleft'], nu2=asfit_raw['kright'], mu=asfit_raw['loc'], sigma=asfit_raw['scale']) for z in x_list]        
        ycdf_shocked= [asymt_cdf(z, alpha=asfit_shocked['skew'], nu1=asfit_shocked['kleft'], nu2=asfit_shocked['kright'], mu=asfit_shocked['loc'], sigma=asfit_shocked['scale']) for z in x_list]
        yzero=asymt_cdf(0, alpha=asfit_raw['skew'], nu1=asfit_raw['kleft'], nu2=asfit_raw['kright'], mu=asfit_raw['loc'], sigma=asfit_raw['scale'])
        yzero_shocked=asymt_cdf(0, alpha=asfit_shocked['skew'], nu1=asfit_shocked['kleft'], nu2=asfit_shocked['kright'], mu=asfit_shocked['loc'], sigma=asfit_shocked['scale'])
        tmp_dic={'AsymT_PDF_x':x_list,'AsymT_PDF_y_before':yvals_raw,'AsymT_CDF_y_before':ycdf_raw,'AsymT_PDF_y_after':yvals_shocked,'AsymT_CDF_y_after':ycdf_shocked}
        dfpdf=pd.DataFrame(tmp_dic)[['AsymT_PDF_x','AsymT_PDF_y_before','AsymT_CDF_y_before','AsymT_PDF_y_after','AsymT_CDF_y_after']]
    #tmp_dic={'Tskew_PDF_x':x_list,'Tskew_PDF_y':yvals,'Tskew_CDF':ycdf}
    #dfpdf=pd.DataFrame(tmp_dic)
        for i,y in enumerate(ycdf_raw):
            if y>0.05:
                q5loc_raw=i
                break
        for i,y in enumerate(ycdf_raw):
            if y>0.1:
                q10loc_raw=i
                break
        
        for i,y in enumerate(ycdf_shocked):
            if y>0.05:
                q5loc_shocked=i
                break
        for i,y in enumerate(ycdf_shocked):
            if y>0.1:
                q10loc_shocked=i
                break     
    
#    for i,x in enumerate(x_list):
#        if x>5:
#            v5loc=i
#            break
        
        titlestr = "Scenario test for "+fitdate.strftime('%m/%d/%Y')+" "+"growth rate"+" forward "+str(horizon)
        lablestr_raw = "Density before shock"
        lablestr_shocked = "Density after shock"
        fig, ax = plt.subplots(1, 1, figsize=(20,10))

#    if fitparam['mode']['constraint']=='Free':
#        loc=tsfit['loc']/tsfit['scale']
#    modx=loc
        ax.set_title(titlestr,fontsize=24)

        if yvals_raw[q5loc_raw]>yvals_shocked[q5loc_raw]:
            
            ax.fill_between(x_list[:q5loc_raw], 0, yvals_raw[:q5loc_raw],  facecolor='c', interpolate=True)
            ax.fill_between(x_list[:q5loc_shocked], 0, yvals_shocked[:q5loc_shocked],  facecolor='r', interpolate=True)
        
        else:
            ax.fill_between(x_list[:q5loc_shocked], 0, yvals_shocked[:q5loc_shocked],  facecolor='r', interpolate=True)
            ax.fill_between(x_list[:q5loc_raw], 0, yvals_raw[:q5loc_raw],  facecolor='c', interpolate=True)
            
        ax.plot(x_list,yvals_raw,'c-',label=lablestr_raw)
        ax.plot(x_list,yvals_shocked,'r-',label=lablestr_shocked)

    
#    if fitparam['mode']['constraint']=='Free':
#        loc=tsfit['loc']/tsfit['scale']
#    modx=loc
#    mody=tskew_pdf(loc, df=tsfit['df'], loc=tsfit['loc'], scale=tsfit['scale'], skew=tsfit['skew'])
#    ax.plot([modx,modx],[0,mody],'r-.')
#    ax.annotate('Mode', xy=(modx, mody),xycoords='data',
#                xytext=(modx-1.5, mody+0.1), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                connectionstyle="arc3"),)
#    ax.plot([x_list[q5loc_raw-1],x_list[q5loc_raw-1]],[0,yq5_raw],'k--')
#    ax.plot([x_list[q5loc_shocked-1],x_list[q5loc_shocked-1]],[0,yq5_shocked],'k--')
#    ax.annotate('GaR 5%, Before', xy=(x_list[max(0,q5loc_raw-1)], yq5_raw), xycoords='data',
#                xytext=(x_list[max(0,q5loc_raw-30)]-2, yq5_raw+0.1), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                                connectionstyle="angle3,angleA=90,angleB=0",facecolor="cyan"),)
#    
#    ax.annotate('GaR 5%, After', xy=(x_list[max(0,q5loc_shocked-1)], yq5_shocked), xycoords='data',
#                xytext=(x_list[max(0,q5loc_shocked-30)]-2, yq5_raw+0.13), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                                connectionstyle="angle3,angleA=90,angleB=0",facecolor="red"),)
#    ax.plot([x_list[q10loc-1],x_list[q10loc-1]],[0,yq10],'k--')
#    ax.annotate('GaR 10%', xy=(x_list[q10loc-1], yq10), xycoords='data',
#                xytext=(x_list[q10loc-1]-1.5, yq10+0.1), textcoords='data',
#                arrowprops=dict(arrowstyle="->",
#                connectionstyle="angle3,angleA=0,angleB=90"),)
        ax.legend(fontsize=24)
        ax.tick_params(labelsize=24)
        plt.ylim(0, max(max(yvals_raw),max(yvals_shocked))+0.2)
        plt.ylabel('Probability Density', fontsize=24)
        plt.xlabel(freq+' GDP(compound annual growth rate)', fontsize=24)    
    
    
        res=[]
        res.append([' ','Before shock','After shock'])
        res.append(['Date of input',fitdate,fitdate])
        res.append(['Horizon forward',horizon,horizon])
        res.append(['Conditional mode',loc_raw,loc_shocked])
        res.append(['GaR5%',x_list[q5loc_raw-1],x_list[q5loc_shocked-1]])
        res.append(['GaR10%',x_list[q10loc_raw-1],x_list[q10loc_shocked-1]])
        res.append(['Growth below 0 probablity',float("{:.4f}".format(yzero)),float("{:.4f}".format(yzero_shocked))])
        res.append(['Skew parameter',asfit_raw['skew'],asfit_shocked['skew']])
        res.append(['Scale',asfit_raw['scale'],asfit_shocked['scale']])   
        return fig,res,dfpdf