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

###############################################################################
#%% Plotting
###############################################################################
## Style of the charts
plt.style.use('seaborn-white')

## Charting parameters : size
from pylab import rcParams


###############################################################################
#%% Partition plot
###############################################################################
def partition_plot(ddatac,dload,group_list,depvar):
    plt.close('all')

    
    ## Prepare the frame of original data
    ddatac = ddatac.set_index(ddatac.date) 
    group_list.sort()
    cr=len(group_list)
    fig, axes = plt.subplots(nrows= 2, ncols=len(group_list), figsize=(11*cr,30))
    
    if len(group_list)>1:
        for g, group in enumerate(group_list):

            group_label = group
        
        ## Upper plot
            ddatac.loc[:,group].plot(ax=axes[0,g])
            axes[0,g].axhline(y=0, c='black', linewidth=0.7)
            axes[0,g].set_title('{} over time'.format(group_label),
                            fontsize=30, y=1.05)
            axes[0,g].set_xlabel('')
            axes[0,g].tick_params(labelsize=18)
            plt.setp(axes[0,g].xaxis.get_majorticklabels(), rotation=70 )
              
        ## Middle plot
            dl1 = dload[dload.group==group]
            dl1['loadings'] = dl1['loadings'] # Now remove absolute
            dl1 = dl1.sort_values(by=['loadings'], ascending=[0])
            sum_abs = np.sum(np.absolute(dl1.loadings))
            dl1['norm_loadings'] = dl1['loadings']/sum_abs
            x_arr = dl1.norm_loadings.values
            y_arr = np.arange(len(dl1.variable))
            axes[1,g].hlines(y_arr, 0, x_arr, color='red')  # Stems
            axes[1,g].plot(x_arr, y_arr, 'D')  # Stem ends
            axes[1,g].plot([0, 0], [y_arr.min(), y_arr.max()], '--')  # Middle bar
            axes[1,g].set_yticks(range(len(dl1.variable.values)))
            ytick=dl1.variable.values
            for ind in range(len(ytick)):
                if len(ytick[ind])>10:
                    ytick[ind]=ytick[ind][0:7]+'...'
            axes[1,g].set_yticklabels(ytick)
            axes[1,g].set_title('Loadings: {}'.format(group_label), fontsize=30)
            axes[1,g].tick_params(labelsize=25)
        ## Increase size to avoid chart to disappear
            ymin, ymax = axes[1,g].get_ylim()
            axes[1,g].set_ylim([ymin-0.5,ymax+0.5])
    elif len(group_list)==1:
            group=group_list[0]
            group_label = group
        
        ## Upper plot
            ddatac.loc[:,group].plot(ax=axes[0])
            axes[0].axhline(y=0, c='black', linewidth=0.7)
            axes[0].set_title('{} over time'.format(group_label),
                            fontsize=30, y=1.05)
            axes[0].set_xlabel('')
            axes[0].tick_params(labelsize=18)
            plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=70 )
              
        ## Middle plot
            dl1 = dload[dload.group==group]
            dl1['loadings'] = dl1['loadings'] # Now remove absolute
            dl1 = dl1.sort_values(by=['loadings'], ascending=[0])
            sum_abs = np.sum(np.absolute(dl1.loadings))
            dl1['norm_loadings'] = dl1['loadings']/sum_abs
            x_arr = dl1.norm_loadings.values
            y_arr = np.arange(len(dl1.variable))
            axes[1].hlines(y_arr, 0, x_arr, color='red')  # Stems
            axes[1].plot(x_arr, y_arr, 'D')  # Stem ends
            axes[1].plot([0, 0], [y_arr.min(), y_arr.max()], '--')  # Middle bar
            axes[1].set_yticks(range(len(dl1.variable.values)))
            ytick=dl1.variable.values
            for ind in range(len(ytick)):
                if len(ytick[ind])>10:
                    ytick[ind]=ytick[ind][0:7]+'...'
            axes[1].set_yticklabels(ytick)
            axes[1].set_title('Loadings: {}'.format(group_label), fontsize=30)
            axes[1].tick_params(labelsize=25)
        ## Increase size to avoid chart to disappear
            ymin, ymax = axes[1].get_ylim()
            axes[1].set_ylim([ymin-0.5,ymax+0.5])
        
        

        
#        ## Correlation plot
#
#        ddg= ddatac.dropna(subset=[group], axis=0).copy()
#        ddg = ddg.set_index(ddg.date)
#        ddg.loc[:,[depvar, group]].plot(ax=axes[2,g], secondary_y=[group])
#
#        axes[2,g].set_title('Correlation with dependent variable : {0}'.format(group_label),
#                            fontsize=30, y=1.05)
#        axes[2,g].tick_params(labelsize=18)
#        axes[2,g].legend(loc='upper center', prop={'size':15})
#        axes[2,g].set_xlabel('')
#        plt.setp(axes[2,g].xaxis.get_majorticklabels(), rotation=70 )
    
    fig.subplots_adjust(hspace=0.55, wspace=0.55)
    fig.suptitle('Supervised data partitioning normalized loadings',
                 y=0.95,fontsize=45)
    #plt.show()
    
    fig1,ax1 = plt.subplots(1, 1, figsize=(15,30))
    gvar={}
    for group in group_list:
        gvar[group]=dload[dload['group']==group]['variance_ratio'].values[0]
    x_arr=[gvar[g] for g in group_list]
    y_arr=np.arange(len(group_list))
    xs=0
    xm=1.2
    ax1.hlines(y_arr, xs, x_arr, color='red')  # Stems
    #ax1.plot([1, 1], [y_arr.min(), y_arr.max()], '--')
    ax1.set_yticks(range(len(group_list)))
    ytick=group_list
    ax1.set_yticklabels(ytick)
    ax1.set_xticks(np.arange(xs,xm,0.2))
    ax1.set_title('Variance ratio of the partition',fontsize=36)
    ax1.tick_params(labelsize=36)
    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim([ymin-0.5,ymax+0.5])
    ax1.grid()
    return (fig,fig1)
