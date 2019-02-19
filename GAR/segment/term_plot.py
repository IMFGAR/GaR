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
def termstruct_plot(df_term,regressors, qlist,hlist):
    plt.close('all')    
    print(qlist)
    print(hlist)
#    qlist.sort()
#    
#    for i in range(len(qlist)):
#        if qlist[i]==0.5:
#            ind05=i
#            break
#    qlist.insert(ind05,'mean')

    ## Variables text
    variable_list_coeff = list(regressors)
    variable_list_coeff.sort()
    #variable_list_label = {'autoregressive':'Autoregressive',
#                           'prices':'Price of Risk',
#                           'quantities': 'Leverage',
#                           'foreign': 'External'}
    
    ## Define the grid
    n=len(variable_list_coeff)
    m=len(hlist)
    if m<=4:
        cs=m
    elif m<=6:
        cs=3
    else:
        cs=4
    rs=m//4+1
    
    termfigs=[]
    
    
    for v, variable in enumerate(variable_list_coeff):

            
        fig = plt.figure(figsize=(cs*8,8*rs+2))
        axes=[]
        gs = GridSpec((n+1)//4+1, min(4,m),hspace=0.35)
        
    ## Plots    
        vs=variable.split('_trans_')
        varn=vs[0]
        if vs[1][-4:]!='None':
            varn+='_'+vs[1]
        for i in range(m):
            axes.append(fig.add_subplot(gs[i//4,i%4]))
                        
            dcv = df_term.loc[(df_term.index == variable),:].copy()
            dcv = dcv.reset_index()
            dcv = dcv.set_index(dcv['quantile'])
            dcv = dcv.reindex(qlist)
            coff='coeff_scale_hz'+str(hlist[i])
            erro='error_hz'+str(hlist[i])
            upper='upper_hz'+str(hlist[i])
            lower='lower_hz'+str(hlist[i])
            erna=dcv[erro].isnull().any()
            # Plot the coefficients
            if erna:
                dcv[erro].plot.bar(color='blue',ax=axes[i])
                x=max(abs(min(dcv[coff].values)),abs(max(dcv[coff].values)))
            else:
                dcv[coff].plot.bar(color='blue',yerr = dcv[erro],ax=axes[i])    
                
                x=max(abs(min(dcv[lower].values)),abs(max(dcv[upper].values)))
                
            axes[i].axhline(y=0, c='black', linewidth=0.7)
            axes[i].set_title('{0}'.format('Horizon '+str(hlist[i])), fontsize=25, y=1.02)
            axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axes[i].set_xlabel('')
            
            axes[i].set_ylim(-x-0.1,x+0.1)
            axes[i].tick_params(labelsize=25)
        
        
        fig.suptitle('Term structure for '+varn, y=1,fontsize=30)
        termfigs.append(fig)
    #plt.text(-19, 1.40, vars_text, fontsize=22, ha='center')
    # fig.savefig('qfit.png')
    # plt.show()

# Plot the R2
    '''
    dcv['R2_in_sample'].plot.bar(color='blue', ax=axes[n])
    axes[n].set_title('R2: ', fontsize=25, y=1.05)
    axes[n].set_xlabel('')
    axes[n].tick_params(labelsize=25)
    axes[n].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axes[n].set_ylim(0,1)
    fig.subplots_adjust(hspace=0.55, wspace=0.25)   
    '''
    return termfigs
