# -*- coding: utf-8 -*-
"""
Tskew distributions and main moments
Based on Zhu and Galbraith JoE 2010 
rlafarguette@imf.org
Time-stamp: "2018-10-05 12:09:12 RLafarguette"
"""


###############################################################################
#%% Modules
###############################################################################
## Numpy
import numpy as np                                      ## Numeric tools
import scipy                                            ## Scientific tools

from scipy.stats import t                               ## Student distribution

## Special functions
from scipy.special import gamma, digamma, polygamma     ## Gamma functions

###############################################################################
#%% Ancillary functions, cf. Zhu and Galbraith JoE 2010
###############################################################################
def K_plain(nu):
    top = gamma((nu+1)/2)
    bottom = np.sqrt(scipy.pi*nu)*gamma(nu/2)
    return(top/bottom)

def alpha_star_plain(alpha, nu1, nu2):
    top = alpha*K(nu1)
    bottom = alpha*K(nu1) + (1-alpha)*K(nu2)
    return(top/bottom)

## To improve speed, vectorize the ancillary functions (used everywhere else)
K = np.vectorize(K_plain, otypes=[np.float64], cache=False)
alpha_star = np.vectorize(alpha_star_plain, otypes=[np.float64], cache=False)

###############################################################################
#%% Expectation of the Assymetric student t, cf. Zhu and Galbraith JoE 2010
###############################################################################
def asymt_mean(alpha=0.5, nu1=1, nu2=1, mu=0, sigma=1):
    astar=alpha_star_plain(alpha,nu1,nu2)
    knu1=K_plain(nu1)
    knu2=K_plain(nu2)
    #Expectation of standard AST
    East=4*(-alpha*astar*nu1*knu1/(nu1-1)+(1-alpha)*(1-astar)*nu2*knu2/(nu2-1))    
    #Scaled by sigma and shift by mu
    ans=East*sigma+mu
    return ans
    
    
###############################################################################
#%% PDF of the Assymetric student t, cf. Zhu and Galbraith JoE 2010
###############################################################################
   
    
def asymt_pdf(y, alpha=0.5, nu1=1, nu2=1, mu=0, sigma=1):
    
    """ 
    Following Zhu and Galbraith, pp. 298 bottom right
    Alpha is the skewness, nu1 and nu2 are the left and right kurtosis
    mu is location (mode) and sigma the scale (variance)
    """
    
    if y <= mu: ## Specify the density on the left tail
        core = (y-mu)/(2*alpha*sigma*K(nu1))
        core2 = np.power(core,2)
        bracket = 1 + (1/nu1)*core2
        bracket_power = np.power(bracket, -(nu1+1)/2)
        pdf = (1/sigma)*bracket_power
        
    else: ## Specify the density on the right tail
        core = (y-mu)/(2*(1-alpha)*sigma*K(nu2))
        core2 = np.power(core, 2)
        bracket = 1 + (1/nu2)*core2
        bracket_power = np.power(bracket, -(nu2+1)/2)
        pdf = (1/sigma)*bracket_power

    return(pdf)

###############################################################################
#%% CDF of the Assymetric student t, cf. Zhu and Galbraith JoE 2010
###############################################################################
def asymt_cdf(y_0, alpha=0.5, nu1=1, nu2=1, mu=0, sigma=1):
    
    """ 
    Following Zhu and Galbraith, pp. 299 top left  
    Alpha is the skewness, nu1 and nu2 are the left and right kurtosis
    mu is location (mode) and sigma the scale (variance)
    """
    
    ## Need to normalize so that it works with unscaled version of cdf
    y = (y_0 - mu)/sigma
    
    left_bracket = np.min([y,0])/(2*alpha_star(alpha, nu1, nu2))
    left_block = 2*alpha*t.cdf(left_bracket, df=nu1)

    right_bracket = np.max([y,0])/(2*(1-alpha_star(alpha, nu1, nu2)))
    right_block = 2*(1-alpha)*(t.cdf(right_bracket, df=nu2) - (1/2))

    cdf = left_block + right_block
    
    return(cdf)

###############################################################################
#%% PPF of the Assymetric student T, cf. Zhu and Galbraith JoE 2010
###############################################################################
def asymt_ppf(p, alpha=0.5, nu1=1, nu2=1, mu=0, sigma=1):

    """ 
    Following Zhu and Galbraith, pp. 299-300  
    Alpha is the skewness, nu1 and nu2 are the left and right kurtosis
    mu is location (mode) and sigma the scale (variance)
    """
       
    left_bracket = np.min([p, alpha])/(2*alpha)
    left_block = 2*alpha_star(alpha, nu1, nu2)*t.ppf(left_bracket, df=nu1)

    right_bracket = (np.max([p,alpha]) + 1 - (2*alpha))/(2*(1-alpha))
    right_block = 2*(1-alpha_star(alpha, nu1, nu2))*(t.ppf(right_bracket, df=nu2))

    ## Need to normalize to get back to an unscaled distribution
    ppf = (left_block + right_block)*sigma + mu
    
    return(ppf)

###############################################################################
#%% Log likelihood of the assymetric student t, cf. Zhu and Galbraith JoE 2010
###############################################################################
def asymt_log_likelihood(sample, theta): # Theta is a coefficient vector
    """ Following Zhu and Galbraith, pp. 300 bottom right """

    ## Unpack
    alpha = theta[0] 
    nu1 = theta[1] 
    nu2 = theta[2] 
    mu = theta[3] 
    sigma = theta[4] 

    def lkl(y, alpha, nu1, nu2, mu, sigma):
        if y <= mu:
            core = (y-mu)/(2*alpha*sigma*K(nu1))
            core2 = np.power(core,2)
            bracket = 1 + (1/nu1)*core2
            log_bracket = np.log(bracket)
            final_block = -((nu1+1)/2)*log_bracket
        else:
            core = (y-mu)/(2*(1-alpha)*sigma*K(nu2))
            core2 = np.power(core,2)
            bracket = 1 + (1/nu2)*core2
            log_bracket = np.log(bracket)
            final_block = -((nu2+1)/2)*log_bracket

        return(final_block)    

    ## To improve speed, vectorize the likelihood
    lkl_vect = np.vectorize(lkl, otypes=[np.float], cache=False)

    ## Apply the vectorize function to the numpy array
    sample_trans = lkl_vect(sample, alpha, nu1, nu2, mu, sigma)
    sum_loglkl = np.sum(sample_trans) # Sum is already vectorized on numpy
    log_lkl = -len(sample)*np.log(sigma) + sum_loglkl
    
    return(log_lkl)

###############################################################################
#%% Hessian Matrix of the log likelihood: Zhu and Galbraith, p 301, on the left
## H(theta) = - I(theta), where I is the Fisher Information Matrix
###############################################################################

## Ancillary functions
def D_plain(nu):
    return(digamma((nu + 1)/2) - digamma((nu/2)))

def D_prime_plain(nu):
    pblock = polygamma(1, ((nu + 1)/2)) - polygamma(1, ((nu/2)))
    return((1/2)*pblock) # Pay attention to derivation rules...

## To improve speed, vectorize the ancillary functions (used everywhere else)
D = np.vectorize(D_plain, otypes=[np.float64], cache=False)
D_prime = np.vectorize(D_prime_plain, otypes=[np.float64], cache=False)

## Hessian is symmetric, just need to get one half of it (still 15 elements...)
def phi_11(alpha, nu1, nu2, mu, sigma):
    block = (nu1 + 1)/(alpha*(nu1 + 3)) + (nu2 + 1)/((1-alpha)*(nu2 + 3))
    return(3*block)

def phi_12(alpha, nu1, nu2, mu, sigma):
    left = (-1)/(nu1 + 1)
    right = (nu1*D(nu1))/(nu1 + 3)
    return(left + right)

def phi_13(alpha, nu1, nu2, mu, sigma):
    left = (1)/(nu2 + 1)
    right = (nu2*D(nu2))/(nu2 + 3)
    return(left - right)

def phi_14(alpha, nu1, nu2, mu, sigma):
    const = (-2)/(3*sigma)
    final = const*phi_11(alpha, nu1, nu2, mu, sigma)
    return(final)

def phi_15(alpha, nu1, nu2, mu, sigma):
    const = 2/sigma
    block = ((nu1)/(nu1 + 3)) - (nu2)/(nu2 + 3)
    final = const*block
    return(final)

def phi_22(alpha, nu1, nu2, mu, sigma):
    const = alpha/2
    block_left = (nu1*np.power(D(nu1),2))/(nu1 + 3)
    block_middle = (2*D(nu1))/(nu1 + 1)
    block_right = D_prime(nu1)
    final = const*(block_left - block_middle - block_right)
    return(final)

def phi_23(alpha, nu1, nu2, mu, sigma):
    return(0)

def phi_25(alpha, nu1, nu2, mu, sigma): 
    final = (alpha/sigma)*phi_12(alpha, nu1, nu2, mu, sigma)
    return(final)

def phi_24(alpha, nu1, nu2, mu, sigma):
    const = 1/sigma
    left_block = 1/(nu1 + 1)
    right_block = ((nu1 + 1)/(nu1 + 3))*D(nu1)
    final = const*(left_block - right_block)
    return(final)

def phi_34(alpha, nu1, nu2, mu, sigma): 
    const = (-1)/sigma
    left_block = 1/(nu2 + 1)
    right_block = ((nu2 + 1)/(nu2 + 3))*D(nu2)
    final = const*(left_block - right_block)
    return(final)

def phi_33(alpha, nu1, nu2, mu, sigma):
    const = (1-alpha)/2
    left_block = (nu2*np.power(D(nu2),2))/(nu2 + 3)
    middle_block = (2*D(nu2))/(nu2+1)
    right_block = D_prime(nu2)
    final = const*(left_block - middle_block - right_block)
    return(final)

def phi_35(alpha, nu1, nu2, mu, sigma):
    const = (alpha-1)/sigma
    final = const*phi_13(alpha, nu1, nu2, mu, sigma)
    return(final)

def phi_44(alpha, nu1, nu2, mu, sigma):
    const = 1/(4*np.power(sigma,2))
    left_0 = (nu1 + 1)/(alpha*(nu1+3))
    left_1 = 1/(np.power(K(nu1),2))
    right_0 = (nu2 + 1)/((1-alpha)*(nu2 + 3))
    right_1 = 1/(np.power(K(nu2),2))
    final = const*((left_0*left_1) + (right_0*right_1))
    return(final)

def phi_45(alpha, nu1, nu2, mu, sigma):
    const = (-2)/(3*sigma)
    final = const*phi_15(alpha, nu1, nu2, mu, sigma)
    return(final)

def phi_55(alpha, nu1, nu2, mu, sigma):
    const = 2/(np.power(sigma,2))
    left = (alpha*nu1)/(nu1 + 3)
    right = ((1-alpha)*nu2)/(nu2 + 3)
    final = const*(left + right)
    return(final)


## Now package everything into a square 5x5 matrix
def asymt_ll_hessian(theta): # Theta is the vector of parameters
    """ This function computes the Hessian matrix of the log_likelihood """
    
    ## Unpack
    alpha = theta[0] 
    nu1 = theta[1] 
    nu2 = theta[2] 
    mu = theta[3] 
    sigma = theta[4] 
    
    ## Create an empty 5x5 matrix
    shape = (5, 5)
    hess = np.empty(shape)
    
    ## Fill it in a symmetric way using the definitions above

    # Diagonal terms
    hess[0,0] = phi_11(alpha, nu1, nu2, mu, sigma)
    hess[1,1] = phi_22(alpha, nu1, nu2, mu, sigma)
    hess[2,2] = phi_33(alpha, nu1, nu2, mu, sigma)
    hess[3,3] = phi_44(alpha, nu1, nu2, mu, sigma)
    hess[4,4] = phi_55(alpha, nu1, nu2, mu, sigma)

    # Off diagonal: "1" family
    hess[0,1] = phi_12(alpha, nu1, nu2, mu, sigma)
    hess[1,0] = phi_12(alpha, nu1, nu2, mu, sigma)

    hess[0,2] = phi_13(alpha, nu1, nu2, mu, sigma)
    hess[2,0] = phi_13(alpha, nu1, nu2, mu, sigma)

    hess[0,3] = phi_14(alpha, nu1, nu2, mu, sigma)
    hess[3,0] = phi_14(alpha, nu1, nu2, mu, sigma)

    hess[0,4] = phi_15(alpha, nu1, nu2, mu, sigma)
    hess[4,0] = phi_15(alpha, nu1, nu2, mu, sigma)

    # Off diagonal: "2" family
    hess[1,2] = phi_23(alpha, nu1, nu2, mu, sigma)
    hess[2,1] = phi_23(alpha, nu1, nu2, mu, sigma)

    hess[1,3] = phi_24(alpha, nu1, nu2, mu, sigma)
    hess[3,1] = phi_24(alpha, nu1, nu2, mu, sigma)

    hess[1,4] = phi_25(alpha, nu1, nu2, mu, sigma)
    hess[4,1] = phi_25(alpha, nu1, nu2, mu, sigma)

    # Off diagonal: "3" family
    hess[2,3] = phi_34(alpha, nu1, nu2, mu, sigma)
    hess[3,2] = phi_34(alpha, nu1, nu2, mu, sigma)

    hess[2,4] = phi_35(alpha, nu1, nu2, mu, sigma)
    hess[4,2] = phi_35(alpha, nu1, nu2, mu, sigma)

    # Off diagonal: "4" family
    hess[3,4] = phi_45(alpha, nu1, nu2, mu, sigma)
    hess[4,3] = phi_45(alpha, nu1, nu2, mu, sigma)

    ## Pay attention that the Hessian is the opposite of the Fisher !!
    return((-1)*hess)
    

