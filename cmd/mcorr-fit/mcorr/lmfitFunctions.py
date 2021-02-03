"""defining lmfit functions for fitCorr.py script"""
from lmfit import Minimizer, Parameters

def d_i(a, theta_i):
    "eq 6"
    d = theta_i/(1+theta_i*a)
    return d

def c_s0(w, a, theta_s, phi_s, f, l):
    "eq 22"
    c_s0 = (1+2*theta_s*a)/(1+2*theta_s*a+phi_s*w*(f+l))
    return c_s0

def c_s1(w, a, phi_s, l, theta_s, f):
    "eq 23"
    c_s1 = (2*phi_s*w*l)/(1+2*theta_s*a+phi_s*w*(f+l))
    return c_s1

def c_s2(phi_s, w, f, l, theta_s, a):
    "eq 24"
    c_s2 = (phi_s*w*(f-l))/(1+2*theta_s*a+phi_s*w*(f+l))
    return c_s2

def Q_p(theta_p, a, phi_p, w, l):
    "eq 25"
    Q_p = 2*((theta_p/(1+theta_p*a))**2)*((1+theta_p*a+phi_p*w*l)/(1+2*theta_p*a+2*phi_p*w*l))
    return Q_p

def residual(pars, x, data=None):
    ##defines the function to be minimized -- the residuals of equation 18 for P2
    ##load in parameters from lmfit
    phi_s = pars["phi_s"]
    f = pars["f"]
    theta_s = pars["theta_s"]
    theta_p = pars["theta_p"]
    w = pars["w"]
    a = pars["a"]
    c_s = pars["c_s"]
    d_theta_p = pars["d_theta_p"]
    d_theta_s = pars["d_theta_s"]
    phi_p = pars["phi_p"]
    ds = pars["d_s"]
    ##define equations to be plugged into eq 18
    cs0 = c_s0(w, a, theta_s, phi_s, f, x)  # eq 22
    d2thetas = d_i(a, 2*theta_s) #eq 20 for 2*theta_s
    cs1 = c_s1(w, a, phi_s, x, theta_s, f)  # eq 23
    dp = d_i(a, theta_p)  # eq 20 for theta_p
    cs2 = c_s2(phi_s, w, f, x, theta_s, a)  # eq 24
    Qp = Q_p(theta_p, a, phi_p, w, x)  # eq 25
    ##FINALLY EQ 18
    Qs = cs0*d2thetas*ds+cs1*ds*dp+cs2*Qp
    P2 = Qs/ds
    if data is None:
        return P2
    return P2 - data

def perform_lmfit(x, y, d_sample):
    """perform least-squares minimization"""
    pfit = Parameters()
    pfit.add(name="phi_s", vary=True, min=0, max=1, value=1e-5)
    pfit.add(name="f", vary=True, min=3, max=3e5, value=1e3)
    pfit.add(name="theta_s", vary=True, min=0, max=d_sample, value=1e-5)
    #pfit.add(name="theta_p", vary=True, min=0, value = 1e-4)
    ##define the fixed params
    pfit.add(name="w", value=2.0/3.0, vary=False)
    pfit.add(name="a", value=4.0/3.0, vary=False)
    pfit.add(name="d_s", value=d_sample, vary=False)
    ##constrained params
    pfit.add(name="c_s", expr="(phi_s*w*f)/(1+theta_s*a+phi_s*w*f)", min=0, max=1) #eq 21
    #pfit.add(name="d_theta_p", expr="theta_p/(1+theta_p*a)", min=0) #eq 20 for theta_p
    pfit.add(name="d_theta_s", expr="theta_s/(1+theta_s*a)", min=0) #eq 20 for theta_s
    #pfit.add(name="d_s", value=d_sample, expr="c_s*d_theta_p+(1-c_s)*d_theta_s", vary=False) ##CONSTRAINING WITH EQ 19
    #pfit.add(name="c_s", expr="(d_s-d_theta_s)/(d_theta_p-d_theta_s)", min=0, max=1) ##CONSTRAINING WITH EQ 19
    ##constraining with eq 26, but re-written ...
    pfit.add(name="theta_p", expr="((1-c_s)*d_theta_s-d_s)/(a*(d_s-d_theta_s)+c_s*(d_theta_s*a-1))", min=0)
    pfit.add(name="phi_p", expr="(theta_p*phi_s)/theta_s", min=0) #eq. 27
    ##collect d_theta_p for later
    pfit.add(name="d_theta_p", expr="theta_p/(1+theta_p*a)", vary=False) #eq 20 for theta_p
    myfit = Minimizer(residual, pfit,
                      fcn_args=(x,), fcn_kws={'data': y},
                      scale_covar=True)

    result = myfit.leastsq()
    return result