# -*- coding: utf-8 -*-
"""
Tskew moments with  Cython optimization
rlafarguette@imf.org
Time-stamp: "2018-10-04 20:26:07 RLafarguette"
"""

###############################################################################
#%% Modules
###############################################################################
import numpy as np                                      ## Numeric methods

from scipy.stats import t                               ## Student distribution
import math                                             ## Math functions
###############################################################################
#%% Probability density function of a TSkew distribution
###############################################################################
def tskew_pdf(x, df, loc, scale, skew):    
    """
    Density function of the tskew distribution 
    Based on the formula in Giot and Laurent (JAE 2003 pp. 650)
    - x = the value to evaluate
    - df: degrees of freedom (>1)
    - location: mean of the distribution
    - scale: standard deviation of the distribution
    - skew: skewness parameter (>0, if ==1: no skew, <1: left skew, >1 right)
    
    NB: I had to parametrize the formula differently to get consistent results
    
    """
    cons = (2/(skew + (1/skew)))/scale
    norm_x =  (x-loc)/scale
    if x < loc :
        pdf = cons*t.pdf(skew*norm_x, df, loc=0, scale=1) # Symmetric t pdf
    elif x >= loc:
        pdf = cons*t.pdf(norm_x/skew, df, loc=0, scale=1) # Symmetric t pdf
    else:
        raise ValueError('Incorrect parameters')

    return(pdf)
    
###############################################################################
#%% Percentage point function (quantile function) of a TSkew distribution
###############################################################################
def tskew_ppf(tau, df, loc, scale, skew):
    """
    Quantile function of the tskew distribution 
    Based on the formula in Giot and Laurent (JAE 2003 pp. 650)
    - tau = the quantile
    - df: degrees of freedom (>1)
    - location: mean of the distribution
    - scale: standard deviation of the distribution (>0)
    - skew: skewness parameter (>0, if ==1: no skew, <1: left skew, >1 right)
    
    NB: I had to parametrize the formula differently (was wrong in their paper)
    """

    threshold = 1/(1+np.power(skew,2))
    if tau < threshold:
        adj_tau = (tau/2)*(1+np.power(skew,2))
        non_stand_quantile = (1/skew)*t.ppf(adj_tau, df=df, loc=0, scale=1)
    elif tau >= threshold:
        adj_tau = ((1-tau)/2)*(1+(1/np.power(skew,2)))
        non_stand_quantile = -skew*t.ppf(adj_tau, df=df, loc=0, scale=1)
    else:
        raise ValueError('Parameters misspecified')
    
    quantile = loc + (non_stand_quantile*scale) # Pay attention to this one !
    
    return(quantile)


###############################################################################
#%% Cumulative distribution of a TSkew distribution:
###############################################################################
def tskew_cdf(x, df, loc, scale, skew):    
    """
    Density function of the tskew distribution 
    Based on the formula in Giot and Laurent (JAE 2003 pp. 650) 
    and Lambert and Laurent (2002) pp. 10
    - x = real value on the support to evaluate
    - df: degrees of freedom (>1)
    - location: mean of the distribution
    - skew: skewness parameter (>0, if ==1: no skew, <1: left skew, >1 right)
    
    NB: I had to parametrize it differently in order to get consistent results
    
    """
    
    sk2 = np.power(skew, 2); inv_sk2 = 1/sk2
    norm_x1 = (x-loc)/scale
    #norm_x2 = x-(loc/scale)
    if x < loc:
        # t.cdf() is the symmetric t cdf
        cdf = (2/(1+sk2))*t.cdf(skew*norm_x1, df, loc=0, scale=1) 
    elif x >= loc:
        cdf = 1 - (2/(1+inv_sk2))*t.cdf(-norm_x1/skew, df, loc=0, scale=1)
    else:
        raise ValueError('Incorrect parameters')

    return(cdf)

###############################################################################
#%% sampling t-skew via inverse transform sampling
###############################################################################
def tskew_sampling(n,tau, df, loc, scale, skew):
    uni=np.random.random_sample(size=n)
    samples=[tskew_ppf(x,tau, df, loc, scale, skew) for x in uni]
    return samples


###############################################################################
#%% get mean value of t-skew
###############################################################################
def tskew_mean(df, loc, scale, skew):
    """
    Note by C. Wang
    Formula from Equation 5 of
    On Bayesian Modeling of Fat Tails and Skewness
    Carmen Fernandez and Mark F. J. Stee, 1998 JASA
    """
    cons1=skew-1/skew
    Mr=math.gamma((df+1)/2)/(math.sqrt(df*math.pi)*math.gamma(df/2))*2*df*scale
    return cons1*Mr+loc

#%%
#
## Test
#
#df0= 1
#loc0= 4
#scale0= 2
#skew0= 5
#
#
#print(tskew_ppf(0.2,df0,loc0,scale0,skew0))
#print(tskew_cdf(tskew_ppf(0.2,df0,loc0,scale0,skew0),df0,loc0,scale0,skew0))



df0= 2
loc0= 1.8
scale0= 3
skew0=0.8
print(tskew_mean(df0, loc0, scale0, skew0),loc0,scale0, skew0)
#
#plt.plot([tskew_pdf(x, df0, loc0, scale0, skew0) for x in np.linspace(-10,10, 200)])
#plt.show()

