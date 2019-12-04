# -*- coding: utf-8 -*-
"""
Partition the data using unsupervised or supervised dimensionality reduction
rlafarguette@imf.org
Time-stamp: "2017-09-27 10:21:59 rlafarguette"
"""

###############################################################################
#%% Modules
###############################################################################
import pandas as pd                                     ## Dataframes
import numpy as np                                      ## Numeric methods
from datetime import datetime as date

## Dimensionality reduction
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from .PLSDA import PLS_DA
from sklearn.preprocessing import scale


## Function definition (zcore)
def zscore(series):
    z = ((series - series.mean())/series.std(ddof=0))
    return(z)
    
###############################################################################
#%% Data partitioning
###############################################################################
class Partition(object):
    """ 
    Partition dataset using either supervised or unsupervised data reduction

    Inputs:
    - data: the dataset to reduce
    - groups_dict: groups of variables to perform the reduction along
    - reduction: either PCA or LDA. If LDA needs to provide a benchmark (str)
    - benchmark: Name of the variable to supervise the LDA reduction with (str)

    Outputs:
    - loadings: the loadings of the reduction (1 if only one variable)
    - partition: results of the partitioning

    Usage:
    Partition(df, gv.groups_dict, reduction='LDA', benchmark='benchmark')

    """
    __description = "Data partitioning using dimensionality reduction"
    __author = "Romain Lafarguette, IMF/MCM, rlafarguette@imf.org"

    ## Initializer
    def __init__(self, data, groups_dict, reduction='PCA',
                 benchmark=None,PLStarget=None):

        ## Parameters
        self.reduction = reduction
        self.benchmark = benchmark
        self.PLStarget = PLStarget
        ## Clean the dataset according to the type of reduction
        if self.reduction == 'LDA':
            if isinstance(self.benchmark, str):
                dc = data.dropna(subset=[self.benchmark], axis=0, how='any')
                #dc = data.dropna(axis=0, how='any').copy()
                self.data = dc.dropna(axis=1, how='any').fillna(method='ffill').copy()
            else:
                raise ValueError('Need a benchmark with supervised reduction')
        elif self.reduction == 'PCA':
            self.data = data.dropna(axis=0, how='all').dropna(axis=1,how='any').fillna(method='ffill').copy()
            
        elif self.reduction == "PLS":
            self.data = data.dropna(axis=0, how='all').dropna(axis=1,how='any').fillna(method='ffill').copy()
            
        else:
            raise ValueError('Reduction parameter misspecified')

        ## Remove constant columns (create problem in the partitioning)
        self.data = self.data.loc[:, self.data.apply(pd.Series.nunique) != 1]
        
        ## Populate the groups only with the variables available in the frame
        self.var_dict = {k:[x for x in groups_dict[k] if x in self.data.columns]
                         for k in groups_dict.keys()}

        ## Estimate the fit
        if self.reduction == 'PCA':
            self.partition_fit_group, self.loading = self.__partition_fit_PCA()

            ## For consistency with the LDA object, rename some instances
            for group in sorted(list(self.partition_fit_group.keys())):
                setattr(self.partition_fit_group[group], 'fit',
                        self.partition_fit_group[group].fit_transform)
            
        elif self.reduction == 'LDA':
            self.partition_fit_group, self.loading = self.__partition_fit_LDA()
            
        elif self.reduction == "PLS":
           self.partition_fit_group, self.loading = self.__partition_fit_PLS()

        else:
            raise ValueError('Reduction parameter misspecified')
            


        # By default, using the original data (can be customized)
        self.partition = zscore(self.partition_data(self.data)) 

    ## Methods
    def __partition_fit_PCA(self):
        """ Run the data partitioning using Principal Component Analysis """
        groups = sorted(list(self.var_dict.keys()))
        pca_fit_group = dict()
        loadings_frame = list()
        
        for group in groups:
            var_list = self.var_dict[group]
            if len(var_list) > 1: # Run the partition
                # Partitionning
                dg = self.data.loc[:, var_list].copy()
                X = scale(dg) # Need to scale the variables before partitioning

                ## Fit the PCA
                pca_fit = PCA(n_components=1).fit(X)
                pca_fit_group[group] = pca_fit

                ## Store the loadings
                dl = pd.DataFrame(pca_fit.components_, index=['loadings'],
                                  columns=var_list).transpose()
                dl['group'] = group
                dl['variance_ratio']=pca_fit.explained_variance_ratio_[0]
                dl['variable'] = dl.index
                loadings_frame.append(dl)
                
            elif len(var_list) == 1: # Loadings are 1
                dl = pd.DataFrame(index=var_list)
                dl['loadings'] = 1
                dl['variance_ratio']=1
                dl['group'] = group
                dl['variable'] = var_list[0]
                loadings_frame.append(dl)

            else: # Empty group: no loading
                dl = pd.DataFrame(columns=['loadings', 'group', 'variable'])
                dl['loadings'] = np.nan
                dl['variance_ratio']=np.nan
                dl['group'] = group
                dl['variable'] = np.nan
                loadings_frame.append(dl)

        dloading = pd.concat(loadings_frame)

        # Return the fit method and the associated loadings                
        return((pca_fit_group, dloading))                        

    def __partition_fit_LDA(self):
        """ Run the data partitioning using Linear Discriminant Analysis """
        groups = sorted(list(self.var_dict.keys()))
        lda_fit_group = dict()
        loadings_frame = list()
        
        for group in groups:
            var_list = self.var_dict[group]
            if len(var_list) > 1: # Run the partition
                # Partitionning
                dg = self.data.loc[:, var_list].copy()
                
                X = scale(dg) # Need to scale the variables before partitioning
                y = self.data.loc[:, self.benchmark].values
                
                ## Fit the LDA using the benchmark
                lda_fit = LDA(n_components=1).fit(X, y)
                lda_fit_group[group] = lda_fit

                ## Store the loadings
                dl = pd.DataFrame(lda_fit.coef_, index=['loadings'],
                                  columns=var_list).transpose()
        
                dl['variance_ratio']=lda_fit.explained_variance_ratio_[0]
                dl['group'] = group
                dl['variable'] = dl.index
                loadings_frame.append(dl)

            elif len(var_list) == 1: # Loadings are 1
                dl = pd.DataFrame(index=var_list)
                dl['loadings'] = 1
                dl['variance_ratio']=1
                dl['group'] = group
                dl['variable'] = var_list[0]
                loadings_frame.append(dl)

            else: # Empty group: no loading
                dl = pd.DataFrame(columns=['loadings', 'group', 'variable'])
                dl['loadings'] = np.nan
                dl['variance_ratio']=np.nan
                dl['group'] = group
                dl['variable'] = np.nan
                loadings_frame.append(dl)

        dloading = pd.concat(loadings_frame)
        
        # Return the fit method and the associated loadings        
        return((lda_fit_group, dloading))        
    
    def __partition_fit_PLS(self):
        """ Run the data partitioning using Principal Component Analysis """
        groups = sorted(list(self.var_dict.keys()))
        pls_fit_group = dict()
        loadings_frame = list()
        
        for group in groups:
            var_list = self.var_dict[group]
            if len(var_list) > 1: # Run the partition
                # Partitionning
                plsdepvar=self.PLStarget[group]
                plsavlreg=var_list
                
                ## Fit the PLS
                pls = PLS_DA(plsdepvar, plsavlreg,self.data)         
                pls_fit = pls.fit
                pls_fit_group[group] = pls_fit

                ## Store the loadings
                dl = pls.summary
                dl['group'] = group
                dl['variable'] = dl.index
                loadings_frame.append(dl)
                
            elif len(var_list) == 1: # Loadings are 1
                dl = pd.DataFrame(index=var_list)
                dl['loadings'] = 1
                dl['vip']=1
                dl['group'] = group
                dl['variable'] = var_list[0]
                loadings_frame.append(dl)

            else: # Empty group: no loading
                dl = pd.DataFrame(columns=['loadings', 'group', 'variable'])
                dl['loadings'] = np.nan
                dl['vip']=np.nan
                dl['group'] = group
                dl['variable'] = np.nan
                loadings_frame.append(dl)

        dloading = pd.concat(loadings_frame)

        # Return the fit method and the associated loadings                
        return((pls_fit_group, dloading))   
    
    def partition_data(self, dataframe):
        """ Return the aggregated data """
        # From the previous step, extract the fitting for each group
        groups = sorted(list(self.var_dict.keys()))
        pfit = self.partition_fit_group # Either PCA or LDA

        ## Prepare to store the data and the loadings
        da = pd.DataFrame(index=dataframe.index)
        
        for group in groups:
            var_list = self.var_dict[group]
            if len(var_list) > 1: # Use the loadings from the partition fit
                dg = dataframe.loc[:, var_list].copy()
                
                # Scale the variables
                X = scale(dg) 
                              
                ## Generate the data using the partitioning fit
                if self.reduction=='PLS':
                    Y = scale(dataframe.loc[:, self.PLStarget[group]].copy())     
                    da[group] = pfit[group].fit_transform(X,Y)[0]                    
                else:
                    da[group] = pfit[group].transform(X)
                                    
                
            elif len(var_list) == 1: # Simply keep the variable as it is
                da[group] = dataframe.loc[:, var_list[0]]
        
            else: # Empty group
                da[group] = np.nan

        return(da)

    
