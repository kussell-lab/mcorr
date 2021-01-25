"""Infer recombination rates by fitting correlation profile"""
from __future__ import print_function
import numpy as numpy
from lmfit import Parameters, Minimizer
from tqdm import tqdm
from . import FitRes

def Power(a, b):
    """compute power"""
    return a**b

def const_r1(x, fBar, phiC, w):
    """calculate r1 assuming constant fragment size"""
    return numpy.where(x < fBar, w*phiC*x, w*phiC*fBar)

def exp_r1(x, fBar, phiC, w):
    """calculate r1 assuming exponetional decay of fragment size"""
    return w*phiC*fBar*(1.0 - numpy.exp(-x/fBar))

def geom_r1(x, fBar, phiC, w):
    """calculate r1 assuming geom distribution"""
    prob = 1.0/fBar
    return w*phiC*fBar*(1.0 - numpy.power(1-prob, x))

def calcP2(thetaS, r1, r2, ds, a):
    """
    calcP2 using expression computed using Mathematica CForm
    """
    v = (2*(r2*thetaS + ds*r1*(1 + r1 + r2 + a*thetaS))* \
        (r2*Power(thetaS,2) + Power(ds,2)*(1 + r1 + r2 + a*thetaS)* \
        (2*Power(r1,2) + r2 + 3*r1*r2 + Power(r2,2) + a*(r1 + 2*r2)*thetaS) - \
        ds*thetaS*(2*r2 + Power(r1 + r2,2) + a*(r1 + 3*r2)*thetaS)))/ \
        (Power(r1 + r2,2)*(1 + 2*r1 + r2 + 2*a*thetaS)* \
        (-(thetaS*(r1 - r2 + a*thetaS)) + ds*(2*r1 + a*thetaS)* \
        (1 + r1 + r2 + a*thetaS)))
    return v

def fcn2min(params, xvalues, yvalues, r1_func):
    """function 2 min"""
    thetaS = params['thetaS']
    phiS = params['phiS']
    f = params['f']
    w = params['w']
    r1 = r1_func(xvalues, f, phiS, w)
    r2 = phiS * w * f - r1
    ds = params['ds']
    a = params['a']
    p2 = calcP2(thetaS, r1, r2, ds, a) / ds
    return p2 - yvalues

def fit_model(xvalues, yvalues, d_sample, r1_func):
    """fitting correlation profile using lmfit"""
    params1 = Parameters()
    params1.add('ds', value=d_sample, vary=False)
    params1.add('thetaS', value=0.00001, min=0, max=d_sample)
    params1.add('f', value=1000, min=3, max=300000)
    ## originally max was 1
    params1.add('phiS', value=0.00005, min=0)
    params1.add('w', value=2.0/3.0, vary=False)
    params1.add('a', value=4.0/3.0, vary=False)
    ##originally thetaP, phiP had no minima
    params1.add('thetaP', expr='(ds*(1 + phiS*w*f + a*thetaS)-thetaS)/ \
                                ((1 - a*ds)*(phiS*w*f + a*thetaS)-(a*ds))', min=0)
    params1.add('phiP', expr='phiS*thetaP/thetaS', min=0)
    params1.add('c', expr='w*phiS*f/(1+w*phiS*f+thetaS*a)')
    params1.add('dp', expr='thetaP/(1+a*thetaP)')
    params1.add('dc', expr='thetaS/(1+a*thetaS)')
    minner1 = Minimizer(fcn2min, params1, fcn_args=(xvalues, yvalues, r1_func))
    try:
        fitres1 = minner1.minimize()
    except:
        fitres1 = None
    return fitres1

def vary_fit(xvalues, yvalues, d_sample, r1_func, f_i, thetaS_i, phiS_i):
    """fitting correlation profile using lmfit"""
    params1 = Parameters()
    params1.add('ds', value=d_sample, vary=False)
    params1.add('thetaS', value=thetaS_i, min=0, max=d_sample)
    params1.add('f', value=f_i, min=3, max=300000)
    ## originally max was 1
    params1.add('phiS', value=phiS_i, min=0)
    params1.add('w', value=2.0/3.0, vary=False)
    params1.add('a', value=4.0/3.0, vary=False)
    ##originally thetaP, phiP had no minima
    params1.add('thetaP', expr='(ds*(1 + phiS*w*f + a*thetaS)-thetaS)/ \
                                ((1 - a*ds)*(phiS*w*f + a*thetaS)-(a*ds))', min=0)
    params1.add('phiP', expr='phiS*thetaP/thetaS', min=0)
    params1.add('c', expr='w*phiS*f/(1+w*phiS*f+thetaS*a)')
    params1.add('dp', expr='thetaP/(1+a*thetaP)')
    params1.add('dc', expr='thetaS/(1+a*thetaS)')
    minner1 = Minimizer(fcn2min, params1, fcn_args=(xvalues, yvalues, r1_func))
    try:
        fitres1 = minner1.minimize()
    except:
        fitres1 = None
    return fitres1

def fit_one(fitdata, r1_func):
    """Fit one data set"""
    xvalues = fitdata.xvalues
    yvalues = fitdata.yvalues
    dsample = fitdata.d_sample
    fitres = fit_model(xvalues, yvalues, dsample, r1_func)
    if fitres is not None:
        try:
            params = fitres.params.valuesdict()
            residual = fitres.residual
        except ZeroDivisionError as error:
            print(error)
            return None
        return FitRes(fitdata.group, residual, params, dsample)
    return None

def fit_p2(fitdatas, r1_func=const_r1, disable_progress_bar=False):
    """Fit p2"""
    all_results = []
    for fitdata in tqdm(fitdatas.getall(), disable=disable_progress_bar):
        fitres = fit_one(fitdata, r1_func)
        if fitres is not None:
            all_results.append(fitres)
    return all_results

