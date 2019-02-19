# -*- coding: utf-8 -*-
"""
Tskew moments with  Cython optimization
rlafarguette@imf.org
Time-stamp: "2018-02-21 16:12:10 RLafarguette"
"""

###############################################################################
#%% Modules
###############################################################################
import pandas as pd
import numpy as np                                      ## Numeric methods
from .tskew import tskew_ppf

from scipy.stats import t                               ## Student distribution
from scipy.stats import norm
from scipy import interpolate
from scipy.optimize import minimize

###############################################################################
#%% T-Skew distance between a set of estimated quantiles and theoretical ones
###############################################################################

def quantile_interpolation(alpha, cond_quant_dict):
    """ 
    Quantile interpolation function, following Schmidt and Zhu (2016) p12
    - Alpha is the quantile that needs to be interpolated
    - cond_quant_dict is the dictionary of quantiles to interpolate on 

    Return:
    - The interpolated quantile
    """

    ## List of quantiles
    qlist = sorted(list(cond_quant_dict.keys()))
    min_q = min(qlist)
    max_q = max(qlist)

    ## Fix the base quantile function (usually a N(0,1))
    base = norm.ppf

    ## Considering multiple cases
    if alpha in qlist: ## No need to interpolate, just on the spot !!
        interp = cond_quant_dict[alpha]

    elif alpha < min_q: ## The left edge
        ## Compute the slope (page 13) 
        b1_up = (cond_quant_dict[max_q] - cond_quant_dict[min_q])
        b1_low = base(max_q) - base(min_q)
        b1 = b1_up/b1_low

        ## Compute the intercept (page 12)
        a1 = cond_quant_dict[min_q] - b1*base(min_q)

        ## Compute the interpolated value
        interp = a1 + b1*base(alpha)
        
    elif alpha > max_q: # The right edge (same formula)
        ## Compute the slope (page 13) 
        b1_up = (cond_quant_dict[max_q] - cond_quant_dict[min_q])
        b1_low = base(max_q) - base(min_q)
        b1 = b1_up/b1_low

        ## Compute the intercept (page 12)
        a1 = cond_quant_dict[min_q] - b1*base(min_q)

        ## Compute the interpolated value
        interp = a1 + b1*base(alpha)

    else: # In the belly
        ## Need to identify the closest quantiles
        local_min_list = [x for x in qlist if x < alpha]
        local_min = max(local_min_list) # The one immediately below

        local_max_list = [x for x in qlist if x > alpha]
        local_max = min(local_max_list) # The one immediately above

        # Compute the slope
        b_up = (cond_quant_dict[local_max] - cond_quant_dict[local_min])
        b_low = base(local_max) - base(local_min)
        b = b_up/b_low

        # Compute the intercept
        a = cond_quant_dict[local_max] - b*base(local_max)
        
        ## Compute the interpolated value
        interp = a + b*base(alpha)

    ## Return the interpolated quantile    
    return(interp)

###############################################################################
#%% Uncrossing 
###############################################################################
def quantile_uncrossing(cond_quant_dict, method='linear'):
    """ 
    Uncross a set of conditional_quantiles using Cherzonukov et al 2010
    Via bootstrapped rearrangement
    
    Input:
    - A dictionary of quantile: conditional quantiles
    - Interpolation method: either linear or probabilistic. 
    The probabilistic quantile interpolation follows Zhu and Schmidt 2016

    Output:
    - A dictionary of quantile: uncrossed conditional quantiles

    """

    ## List of quantiles
    ql = sorted(list(cond_quant_dict.keys()))
    cond_quant = [cond_quant_dict[q] for q in ql] # Because dict is not ordered
    np.random.seed(2018)
    ## Check if the quantiles are crossing in the first place
    if sorted(cond_quant) == cond_quant:
        print('Conditional quantiles already sorted !')
        cond_quant_uncrossed_dict = cond_quant_dict
    else:
        if method=='linear':         
            ## Use a linear interpolation for the quantile function
            inter_lin = interpolate.interp1d(ql, cond_quant,
                                             fill_value='extrapolate')

            ## Bootstrap the quantile function
            bootstrap_qf = inter_lin(np.random.uniform(0,1,1000))

            ## Now compute the percentiles of the bootstrapped quantiles 
            cond_quant_uncrossed = [np.percentile(bootstrap_qf, 100*q)
                                    for q in ql]

            ## They are the uncrossed quantiles !
            cond_quant_uncrossed_dict = dict(zip(ql, cond_quant_uncrossed))

        elif method=='probabilistic':
            ## Use Schmidt and Zhu (2016) approach
            bootstrap_qf = [quantile_interpolation(u, cond_quant_dict)
                            for u in np.random.uniform(0,1,1000)]

            ## Now compute the percentiles of the bootstrapped quantiles 
            cond_quant_uncrossed = [np.percentile(bootstrap_qf, 100*q)
                                    for q in ql]

            ## They are the uncrossed quantiles !
            cond_quant_uncrossed_dict = dict(zip(ql, cond_quant_uncrossed))

        else:
            raise ValueError('Interpolation method misspecified')
            
    ## Return the uncrossed quantiles    
    return(cond_quant_uncrossed_dict)


def tskew_distance(quantile_list, cond_quant, 
                   df, loc, scale, skew):
    """ Return the distance between theoretical and actual quantiles"""

    def tskew_tau(tau):
        """ Function which only depends on a given tau """
        return(tskew_ppf(tau, df=df, loc=loc, scale=scale, skew=skew))

    tskew_ppf_vectorized = np.vectorize(tskew_tau, otypes=[np.float])
    
    theoretical_quant = tskew_ppf_vectorized(quantile_list)  

    diff = np.subtract(theoretical_quant, cond_quant)
    diff2 = np.power(diff,2)    
    msse = np.sum(diff2)
    
    loc_tskew=loc/scale
    for i in range(len(quantile_list)):
        if quantile_list[i]==0.25:
            lowq=cond_quant[i]
        if quantile_list[i]==0.75:
            highq=cond_quant[i]
    alpha=10
    if loc_tskew<=highq and loc_tskew>=lowq:
        penalty=0
    else:
        penalty=alpha*min((lowq-loc_tskew)**2,(highq-loc_tskew)**2)
    mssepen=msse+penalty
    return(mssepen)

###############################################################################
#%% Optimal TSkew fit based on a set of conditional quantiles and a location
###############################################################################
def tskew_fit(conditional_quantiles, fitparams):
    """ 
    Optimal TSkew fit based on a set of conditional quantiles and a location
    Inputs:
        - conditional_quantiles (dictionary): quantiles & conditional value
        - loc: location. Can be estimated as a conditional mean via OLS

    Output:
        - A dictionary with optimal scale and skewness, as well as df and loc 
    """
    conditional_quantiles=quantile_uncrossing(conditional_quantiles)
    quantile_list = np.sort(list(conditional_quantiles.keys()))
    
    ######################
    #Generate Parameters##
    ######################
    if fitparams['skew_low']!='Free':
        pass
    

    ## Interquartile range (proxy for volatility)
    try:
        IQR = np.absolute(conditional_quantiles[0.75] - conditional_quantiles[0.25])
        IQR = np.clip(IQR, 1, 10) # Avoid degenerate interquartile range
        # At least 1 pp growth and at most 10 ppt growth in the interquartile
    except:
        raise ValueError('Need to provide estimate for 25% and 75% quantiles')


    ## Upper-bound for the scale
#    scale_up = IQR/1.63 + 0.2# When skew=1, variance exactly = IQR/1.63 
#    scale_down = np.sqrt(IQR)/2 +0.1# Good lower bound approximation
    
    # Avoid crossing: sort
    cond_quant = np.sort(list(conditional_quantiles.values()))
    
    ## Define the boundaries of the conditional quantiles
    #min_o = np.nanmin(cond_quant)
    #max_o = np.nanmax(cond_quant)

    ## Conditional mean can not be inside the conditional quantiles
    ## Else the distribution would be completely degenerated

    # cond_mean = np.clip(loc, min_o, max_o)
    #print(fitparams)
    if fitparams['var_low']['constraint']=='Fixed':
        scale_down=fitparams['var_low']['value']
    else:
        scale_down = np.sqrt(IQR)/2 +0.1# Good lower bound approximation
        
    if fitparams['var_high']['constraint']=='Fixed':
        scale_up=fitparams['var_high']['value']            
    else:
        scale_up = IQR/1.63 + 0.2# When skew=1, variance exactly = IQR/1.63 
        
    if fitparams['skew_low']['constraint']=='Fixed':
        skew_low =fitparams['skew_low']['value']
    else:
        skew_low = 0.1 # Default lower bound approximation
        
    if fitparams['skew_high']['constraint']=='Fixed':
        skew_high=fitparams['skew_high']['value']
        x0_f = [(scale_down+scale_up)/2, (skew_low + skew_high)/2] # Initial values            
    else:
        skew_high = 3 # Default higher bound approximation
        x0_f = [IQR/1.63 + 0.1, 1]
    
    if fitparams['mode']['constraint']=='Fixed':
        loc=fitparams['mode']['value']
    else:
        loc=conditional_quantiles[0.5]
    
    print( scale_down, scale_up, skew_high,skew_low)
    
    
    if fitparams['dof']['constraint']=='Fixed':        
        o_df = fitparams['dof']['value'] # Degrees of freedom
    else:
        o_df = 2


    ## Two values optimizer: on both conditional variance and skewness
    def mult_obj_distance(x): # x is a vector
        """ Multiple parameters estimation """
        ## Unpack the vector
        scale = x[0]
        skew = x[1]

        # Run the optimizer
        obj = tskew_distance(quantile_list=quantile_list,
                             cond_quant=cond_quant,
                             df=o_df, loc=cond_mean, scale=scale, skew=skew)
        return(obj)
    
    
    
    ## Two values optimizer: on both conditional variance and skewness
    def mult_obj_distance3(x): # x is a vector
        """ Multiple parameters estimation """
        ## Unpack the vector
        scale = x[0]
        skew = x[1]
        dloc = x[2]
        # Run the optimizer
        obj = tskew_distance(quantile_list=quantile_list,
                             cond_quant=cond_quant,
                             df=o_df, loc=dloc, scale=scale, skew=skew)
        return(obj)
    
    
    ## Run the optimizer
    locs = loc+0.5
    cond_mean=0
    cdmeanmax=loc+10
    cdmeanmin=max(0.01,loc-10)
    #print(fitparams['mode']['constraint'])
    
    
    if fitparams['mode']['constraint']!='Free':
    #bisection optimize for location
        maxit=0
        while maxit<100 and abs(locs-loc)>0.00001:
            
            cond_mean=(cdmeanmin+cdmeanmax)/2
                
            #print(locs,loc,cdmeanmin,cdmeanmax)        
    
    # Fix the boundaries to avoid degenerative distributions
            bnds_f = ((scale_down, scale_up), (skew_low , skew_high))
            #print(bnds_f)
            res = minimize(mult_obj_distance, x0=x0_f,
                           bounds=bnds_f, method='SLSQP',
                           options={'maxiter':1000,  'ftol': 1e-04, 'eps': 1.5e-06})
        
            o_scale, o_skew  = res.x
            locs=cond_mean/o_scale
            if locs>loc:
                cdmeanmax=cond_mean
            else:
                cdmeanmin=cond_mean
            maxit+=1
            #print(loc,locs,cond_mean)
        
    ## Package the results into a dictionary
        fit_dict = {'loc': float("{:.4f}".format(cond_mean)),
                    'df': int(o_df),
                    'scale': float("{:.4f}".format(o_scale)),
                    'skew': float("{:.4f}".format(o_skew))}
    
        return(fit_dict)
    else:
        x0_f.append(0) # Initial values
    
        # Fix the boundaries to avoid degenerative distributions
        bnds_f = ((scale_down, scale_up),  (skew_low , skew_high), (-20,20))
        #print(bnds_f)
        res = minimize(mult_obj_distance3, x0=x0_f,
                       bounds=bnds_f, method='SLSQP',
                       options={'maxiter':1000,  'ftol': 1e-04, 'eps': 1.5e-06})
        
        o_scale, o_skew, o_loc  = res.x
        
    ## Package the results into a dictionary
        fit_dict = {'loc': float("{:.4f}".format(o_loc)),
                    'df': int(o_df),
                    'scale': float("{:.4f}".format(o_scale)),
                    'skew': float("{:.4f}".format(o_skew))}
        #print(fit_dict)
        return(fit_dict)


###############################################################################
###############################################################################
#%% Risk analysis : Delta VaR under different scenarios
###############################################################################
###############################################################################
def delta_VaR(quant_model, central_scenario, simulated_scenario):
    """ 
    Compute the variation of the VaR from a central scenario 

    Inputs:
    - quant_model: a quantile_reg object, including the regression coefficients
    - central scenario: values of the covariates to estimate the model on
    - simulated scenario: another set of covariates to measure the delta risk
    
    Output:
    - the parameters of the tskew-fit under both scenarios

    """

    ## Shorter names
    qr = quant_model
    
    def opt_parameters(scenario):
        """ Return the optimal parameters associated with a scenario """
        
        ## Compute the associated simulated quantiles
        sq = qr.cond_quantiles(predictors=scenario) 

        ## Retrieve the conditional quantiles, conditional mean and df
        quantile_list = qr.quantile_list
        cond_quant = sq.loc[sq['tau'] != 'mean','conditional_quantile_mean']
        cond_mean = sq.loc[sq['tau'] == 'mean','conditional_quantile_mean'].values[0]
        cq_variables = ['conditional_quantile_mean',
                        'conditional_quantile_mean_ci_lower',
                        'conditional_quantile_mean_ci_upper']
        cond_var = np.nanvar(sq.loc[:, cq_variables].values)
        degree_freedom = qr.data.shape[0] - len(qr.regressors) - 1

        ## Estimate the tskew_fit associated with the conditional quantiles
        tsk = tskew_fit(quantile_list=quantile_list,
                        conditional_quantiles=cond_quant,
                        df=degree_freedom, loc=cond_mean, cond_var=cond_var)

        return(tsk)


    ## Compute the optimal parameters for both scenarios
    c_tsk = opt_parameters(central_scenario) # Central
    s_tsk = opt_parameters(simulated_scenario) ## Simulated

    ## Package them into a frame
    dc = pd.DataFrame.from_dict(c_tsk, orient='index').transpose()
    dc.index = ['central']
    dc.insert(0,'scenario', 'central')

    ds = pd.DataFrame.from_dict(s_tsk, orient='index').transpose()
    ds.index = ['simulated']
    ds.insert(0,'scenario', 'simulated')

    dfinal = pd.concat([dc, ds], axis=0)
    
    return(dfinal)

    
