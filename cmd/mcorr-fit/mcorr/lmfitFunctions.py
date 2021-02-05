"""defining lmfit functions for fitCorr.py script"""
from lmfit import Minimizer, Parameters
import numpy as np

def c_s(phi_s, w, f, theta_s, a):
    """eq 21"""
    c_s = (phi_s*w*f)/(1+theta_s*a+phi_s*w*f)
    return c_s

def d_i(a, theta_i):
    "eq 6"
    d = theta_i/(1+theta_i*a)
    return d

def c_s0(c_s1, c_s2, l):
    "eq 14"
    return np.ones(len(l))-c_s1-c_s2

def c_s1(w, a, phi_s, l, theta_s, f):
    "eq 23"
    ##for l < f
    c_s1less = (2*phi_s*w*l)/(1+2*theta_s*a+phi_s*w*(f+l))
    ## for f >= l
    c_s1greater = (2*phi_s*w*l)/(1+2*theta_s*a+phi_s*w)
    return np.where(l < f, c_s1less, c_s1greater)

def c_s2(phi_s, w, f, l, theta_s, a):
    "eq 24"
    c_s2 = (phi_s*w*(f-l))/(1+2*theta_s*a+phi_s*w*(f+l))
    return np.where(l < f, c_s2, 0)

def Q_p(theta_p, a, phi_p, w, l):
    "eq 25"
    Q_p = 2*((theta_p/(1+theta_p*a))**2)*((1+theta_p*a+phi_p*w*l)/(1+2*theta_p*a+2*phi_p*w*l))
    return Q_p

def residual(pars, x, data=None):
    "defines the function to be minimized -- the residuals of equation 18 for P2"
    ##load in parameters from lmfit
    phi_s = pars["phi_s"]
    #f = pars["f"]
    theta_s = pars["theta_s"]
    theta_p = pars["theta_p"]
    w = pars["w"]
    a = pars["a"]
    c_s = pars["c_s"]
    #d_theta_p = pars["d_theta_p"]
    d_theta_s = pars["d_theta_s"]
    phi_p = pars["phi_p"]
    ds = pars["d_s"]
    f = pars["f"]

    ##define equations to be plugged into eq 18
    #cs0 = c_s0(w, a, theta_s, phi_s, f, x) # eq 22
    d2thetas = d_i(a, 2*theta_s) #eq 20 for 2*theta_s
    #cs1 = pars["c_s1"]
    cs1 = c_s1(w, a, phi_s, x, theta_s, f) # eq 23
    dp = d_i(a, theta_p) # eq 20 for theta_p
    cs2 = c_s2(phi_s, w, f, x, theta_s, a) # eq 24
    cs0 = c_s0(cs1, cs2, x)
    Qp = Q_p(theta_p, a, phi_p, w, x) # eq 25
    ##FINALLY EQ 18
    Qs = cs0*d2thetas*ds+cs1*ds*dp+cs2*Qp
    P2 = Qs/ds
    if data is None:
        return P2
    return P2 - data


def perform_lmfit(x, y, d_sample):
    "perform the fitting with lmfit"
    pfit = Parameters()
    pfit.add(name="phi_s", vary=True, min=0, value=1e-5) ##originally had upper bound of 1
    pfit.add(name="f", vary=True, value=7.5e2, min=3) ##originally min=3; value=1e3
    pfit.add(name="theta_s", vary=True, min=0, value=1e-5)
    ##define the fixed params
    pfit.add(name="w", value=2.0/3.0, vary=False)
    pfit.add(name="a", value=4.0/3.0, vary=False)
    pfit.add(name="d_s", vary=False, value=d_sample)
    ##constrained params
    ##originally 0 to 1 for c_s
    pfit.add(name="c_s", expr="(phi_s*w*f)/(1+theta_s*a+phi_s*w*f)") #eq 21
    pfit.add(name="d_theta_s", expr="theta_s/(1+theta_s*a)") #eq 20 for theta_s (for eq 26)
    pfit.add(name="theta_p", expr="((1-c_s)*d_theta_s-d_s)/(a*(d_s-d_theta_s)+c_s*(d_theta_s*a-1))") #eq 26
    pfit.add(name="phi_p", expr="(theta_p*phi_s)/theta_s") #eq. 27
    pfit.add(name="d_theta_p", expr="theta_p/(1+theta_p*a)") #eq 20 for theta_p (for outputs)

    #do the fitting
    myfit = Minimizer(residual, pfit,
                      fcn_args=(x,), fcn_kws={'data': y},
                      scale_covar=True)
    ##pick least squares with TRF
    result = myfit.least_squares()
    return result