# -*- coding: utf-8 -*-
"""
Useful plotting functions
rlafarguette@imf.org
Time-stamp: "2018-02-05 22:35:15 RLafarguette"
Editted by cwang2@imf.org
"""

###############################################################################
#%% Modules
###############################################################################
#import os, sys, importlib                             ## Operating system
import pandas as pd                                   ## Dataframes
import numpy as np                                    ## Numeric tools
import matplotlib.pyplot as plt                       ## Plotting
import seaborn as sns                                 ## Plotting
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
###############################################################################
#%% Plotting
###############################################################################
## Style of the charts
plt.style.use('seaborn-white')

## Charting parameters : size
from pylab import rcParams
plt.close('all')  

###############################################################################
#%% Coefficients plotting
###############################################################################
def coeff_plot(dcoeffc, regressors, qlist,nhz):
    plt.close('all')    
    qlist.sort()
    
    for i in range(len(qlist)):
        if qlist[i]==0.5:
            ind05=i
            break
    qlist.insert(ind05,'mean')

    ## Variables text
    variable_list_coeff = list(regressors)
    variable_list_coeff.sort()
    #variable_list_label = {'autoregressive':'Autoregressive',
#                           'prices':'Price of Risk',
#                           'quantities': 'Leverage',
#                           'foreign': 'External'}
    
    ## Define the grid
    n=len(variable_list_coeff)

    fig, ax = plt.subplots(n, 1, figsize=(20,9*n))

    ## Plots    
    colorlist=['red','blue','green','cyan','magenta','orange','lime','violet','crimson']
    inds=np.arange(len(qlist))
    bar_width=0.1
    for v, variable in enumerate(variable_list_coeff):
        vs=variable.split('_trans_')
        varn=vs[0]
        if vs[1][-4:]!='None':
            varn+='_'+vs[1]
        if len(varn)>20:
            variable_label = varn[:17]+'...'
        else:
            variable_label = varn
        maxv=-99999999                    
        for hind in range(nhz):
            cn=[]
            for q in qlist:
                cn.append(dcoeffc[(dcoeffc.index==variable) & (dcoeffc['quantile']==q)]['coeff_scale_PROJ'+str(hind+1)].values[0])
            maxv=max(maxv,max([abs(a) for a in cn]))
            if n>1:
                ax[v].bar(inds+hind*bar_width,cn,bar_width,alpha=0.7,color=colorlist[hind],label='PROJ'+str(hind+1))
                plt.sca(ax[v])
            else:
                ax.bar(inds+hind*bar_width,cn,bar_width,alpha=0.7,color=colorlist[hind],label='PROJ'+str(hind+1))
            plt.xticks(inds+nhz//2*bar_width,qlist)
        '''    
        dcv = dcoeffc.loc[(dcoeffc.variable == variable),:].copy()
        dcv = dcv.reset_index()
        dcv = dcv.set_index(dcv['quantile'])
        dcv = dcv.reindex(qlist)
        erna=dcv['errors'].isnull().any()
            # Plot the coefficients
        if erna:
            dcv['coeff_scale'].plot.bar(color='blue',ax=axes[v])
            x=max(abs(min(dcv['coeff_scale'].values)),abs(max(dcv['coeff_scale'].values)))
        else:
            dcv['coeff_scale'].plot.bar(color='blue',yerr = dcv.errors,ax=axes[v])
            x=max(abs(min(dcv['lower'].values)),abs(max(dcv['upper'].values)))
        '''
        if n>1:
            ax[v].axhline(y=0, c='black', linewidth=0.7)
            ax[v].set_title('{0}'.format(variable_label), fontsize=25, y=1.05)
            ax[v].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax[v].set_xlabel('')
            ax[v].tick_params(labelsize=25)  
            ax[v].legend(fontsize=20,bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
            ax[v].set_ylim(-1.1*maxv,1.1*maxv)
        else:
            ax.axhline(y=0, c='black', linewidth=0.7)
            ax.set_title('{0}'.format(variable_label), fontsize=25, y=1.05)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.set_xlabel('')
            ax.tick_params(labelsize=25)  
            ax.legend(fontsize=20,bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
            ax.set_ylim(-1.1*maxv,1.1*maxv)
    if n>1:   
        fig.suptitle('Quantile regressions coefficients',y=0.92,fontsize=30)
    else:
        fig.suptitle('Quantile regressions coefficients',y=1,fontsize=30)
    
    #plt.text(-19, 1.40, vars_text, fontsize=22, ha='center')
    # fig.savefig('qfit.png')
    # plt.show()

    return(fig)
