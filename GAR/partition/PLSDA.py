# -*- coding: utf-8 -*-
"""
Run a PLS - Discriminant Analysis on a set of variables and target variables
Romain Lafarguette, rlafarguette@imf.org
Time-stamp: "2019-03-24 02:06:08 Romain"
"""

###############################################################################
#%% Modules and methods
###############################################################################
## Modules imports
import pandas as pd                                     ## Data management
import numpy as np                                      ## Numeric tools

## Method imports
from sklearn.cross_decomposition import PLSRegression   ## PLS

###############################################################################
#%% Ancillary functions 
###############################################################################
def zscore(series):
    """ Return the Z-score of a series """
    zs = (series-series.mean())/series.std(ddof=0)
    return(zs)


def pls_reduction(depvars, regvars, df):
    assert isinstance(depvars, list), 'Dependent variable(s) should be in list'
    avl_regs = [x for x in regvars if x in df.columns]
    pls_series = PLS_DA(depvars, avl_regs, df).component
    return(pls_series)


def num_days(dates_tuple):
    """ Return the number of days in a tuple """
    min_ = pd.to_datetime(min(dates_tuple))
    max_ = pd.to_datetime(max(dates_tuple))
    return((max_ - min_).days)

###############################################################################
#%% PLS Discriminant Analysis Class Wrapper
###############################################################################
class PLS_DA(object):
    """ 
    Data reduction through PLS-discriminant analysis and variables selection 

    Parameters
    ----------
    dep_vars : list; list of dependent variables
    reg_vars : list; list of regressors variables
    data : pandas df; data to train the model on
    num_vars : 'all', integer; number of variables to keep, ranked by VIP
        if 'all': keep all the variables
    
    Return
    ------
    first_component : the first component of the PLS of the Xs reduction
    output_frame : frame containing the variables and their transformation
    summary_frame : frame with the results of the model (loadings, vip, R2)

    """
    __description = "Partial Least Squares with variables selection"
    __author = "Romain Lafarguette, IMF, rlafarguette@imf.org"

    #### Class Initializer
    def __init__(self, dep_vars, reg_vars, data, num_vars='all'):

        #### Attributes
        self.dep_vars = dep_vars
        self.reg_vars = reg_vars
        self.df = data.dropna(subset=self.reg_vars)

        ## Put parametrized regression as attribute for consistency
        self.pls1 = PLSRegression(n_components=1, scale=True) # Always scale

        ## Unconstrained fit: consider all the variables 
        self.ufit = self.pls1.fit(self.df[self.reg_vars],
                                  self.df[self.dep_vars])

        ## Return the component and summary of the unconstrained model
        ## To save computation time, run it by default for both models        
        self.component_unconstrained = self.__component(self.ufit,
                                                        self.dep_vars,
                                                        self.reg_vars, self.df)

        self.target_unconstrained = self.__target(self.ufit,
                                                  self.dep_vars,
                                                  self.reg_vars, self.df)

        self.summary_unconstrained = self.__summary(self.ufit, self.dep_vars,
                                                    self.reg_vars, self.df)

        ## Variables selection
        if num_vars == 'all': # Unconstrained model: constrained is identical
            self.top_vars = self.reg_vars # The best variables are the full set
            self.fit = self.ufit
            self.component = self.component_unconstrained
            self.target = self.target_unconstrained
            self.summary = self.summary_unconstrained
            
        elif num_vars > 0: ## Constrained model
            self.num_vars = int(num_vars)
            
            ## Identify the most informative variables from the unconstrained
            self.top_vars = list(self.summary_unconstrained.sort_values(
                by=['vip'], ascending=False).index[:self.num_vars])

            ## Now run the constrained fit on these variables
            self.cfit = self.pls1.fit(self.df[self.top_vars],
                                      self.df[self.dep_vars])

            ## Return the main attributes, consistent names with unconstrained
            self.fit = self.cfit
            
            self.component = self.__component(self.cfit, self.dep_vars,
                                              self.top_vars, self.df)
            
            self.target = self.__target(self.cfit, self.dep_vars,
                                        self.top_vars, self.df)

            self.summary = self.__summary(self.cfit, self.dep_vars,
                                          self.top_vars, self.df)
                      
        else:
            raise ValueError('Number of variables parameter misspecified')

        
    #### Internal class methods (start with "__")
    def __vip(self, model):
        """ 
        Return the variable influence in the projection scores
        Input has to be a sklearn fitted model
        Not available by default on sklearn, so it has to be coded by hand
        """
        ## Get the score, weights and loadings
        t = model.x_scores_
        w = model.x_weights_
        q = model.y_loadings_
        p, h = w.shape

        ## Initialize the VIP
        vips = np.zeros((p,))
        s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
        total_s = np.sum(s)

        for i in range(p):
            weight = [(w[i,j] / np.linalg.norm(w[:,j]))**2 for j in range(h)]
            vips[i] = np.sqrt(p*(s.T @ weight)/total_s)
        return(vips)

    def __summary(self, fit, dep_vars, reg_vars, df):
        """
        Return the summary information about the fit
        """
        
        ## Store the information into a pandas dataframe
        dr = pd.DataFrame(reg_vars, columns=['variable'], index=reg_vars)
        dr['loadings'] = fit.x_loadings_ # Loadings
        dr['vip'] = self.__vip(fit) ## Variable importance in projection
        dr['score'] = fit.score(df[reg_vars],df[dep_vars]) # Score
        
        ## Return the sorted summary frame
        return(dr.sort_values(by=['vip'], ascending=False))
    
    ## Write short ancillary functions to export the results into pandas series
    def __component(self, fit, dep_vars, reg_vars, df):
        """
        Return the first component of the fit
        """
        comp = fit.fit_transform(df[reg_vars], df[dep_vars])[0]
        comp_series = pd.Series(comp.flatten(), index=self.df.index)
        return(comp_series)

    def __target(self, fit, dep_vars, reg_vars, df):
        """
        Return the target of the fit (reduced in case of multiple variables)
        """
        target = fit.fit_transform(df[reg_vars], df[dep_vars])[1]
        target_series = pd.Series(target.flatten(), index=self.df.index)
        return(target_series)

    
    #### Standard class methods (no "__")
    def predict(self, dpred):
        """ 
        Apply the dimension reduction learned on new predictors
        Input:
            - dpred: Pandas frame with the predictors 

        Output:
            - Reduced dataframe using the same loadings as estimated in-sample
 
        """
        
        ## Need to select exactly the predictors which have been estimated
        dp = dpred[self.top_vars].dropna()

        ## Run the projection
        dproj = pd.Series(self.fit.predict(dp).flatten(), index=dp.index)

        ## Scaling pb: prediction and fit don't match in sample (they should) !
        # Create the in-sample prediction
        dproj_in = pd.Series(self.fit.predict(self.df[self.top_vars]).flatten(),
                             index=self.df.index)
        
        # Adjust based on the in-sample projection !
        mean_adj =  dproj_in.mean() - self.component.mean()
        scale_adj = dproj_in.std()/self.component.std()

        # Nicely adjusted ! 
        dproj_mod = (dproj-mean_adj)/scale_adj
        
        return(dproj_mod)






