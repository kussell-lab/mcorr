"""Infer recombination rates by fitting correlation profile"""
from __future__ import print_function
import numpy as numpy
from lmfit import Parameters, Minimizer
from tqdm import tqdm
from . import FitRes

def Power(a, b):
    """compute power"""
    return a**b

def const_r1(x, fBar, phiC):
    """calculate r1 assuming constant fragment size"""
    return numpy.where(x < fBar, phiC*x, phiC*fBar)

def exp_r1(x, fBar, phiC):
    """calculate r1 assuming exponetional decay of fragment size"""
    return phiC*fBar*(1.0 - numpy.exp(-x/fBar))

def geom_r1(x, fBar, phiC):
    """calculate r1 assuming geom distribution"""
    prob = 1.0/fBar
    return phiC*fBar*(1.0 - numpy.power(1-prob, x))

def calcP2(fBar, thetaC, phiC, d, x):
    """
    calcP2 using expression from Mathematica
    Yes! The line is super long!
    """
    r1 = const_r1(x, fBar, phiC)
    r2 = phiC * fBar - r1
    v = (4*(2*d*r1*(1 + r1 + r2) + (r1 + 3*d*r1 + 3*r2)*thetaC)*(Power(2*d*(1 + r1 + r2) + 3*d*thetaC,2)*(8*Power(r1,2) + 9*r1*(2*r2 + thetaC) + 9*r2*(1 + r2 + 3*thetaC)) + 2*Power(thetaC,2)*(-2*Power(r1,2) + 9*r1*(2*r2 - thetaC) + 9*r2*(2 + 2*r2 + 3*thetaC)) + d*thetaC*(2*(1 + r1 + r2) + 3*thetaC)*(4*Power(r1,2) - 9*r1*(2*r2 + thetaC) - 9*r2*(4 + 2*r2 + 9*thetaC))))/(Power(r1 + r2,2)*(3 + 4*r1 + 3*r2 + 9*thetaC)*(6 + 4*r1 + 6*r2 + 9*thetaC)*(d*(2*(1 + r1 + r2) + 3*thetaC)*(8*r1 + 9*thetaC) - 2*thetaC*(2*r1 - 6*r2 + 9*thetaC)))
    return v

def fcn2min(params, xvalues, yvalues):
    """function 2 min"""
    fbar = params['fbar']
    dsample = params['dsample']
    phi_clonal = params['phi_clonal']
    theta_clonal = params['theta_clonal']
    p2 = calcP2(fbar, theta_clonal, phi_clonal, dsample, xvalues) / dsample
    return p2 - yvalues

def fit_model(xvalues, yvalues, d_sample):
    """Do fitting using the Model 1"""
    params1 = Parameters()
    params1.add('dsample', value=d_sample, vary=False)
    params1.add('theta_clonal', value=0.00001, min=0, max=d_sample)
    params1.add('fbar', value=1000, min=3, max=300000)
    params1.add('phi_clonal', value=0.0005, min=0, max=1)
    params1.add('theta', expr='(-theta_clonal + dsample*(1 + fbar*phi_clonal + (3*theta_clonal)/2.))/((-3*dsample)/2. + (1 - (3*dsample)/2.)*(fbar*phi_clonal + (3*theta_clonal)/2.))')
    params1.add('phi', expr='(phi_clonal*(-theta_clonal + dsample*(1 + fbar*phi_clonal + (3*theta_clonal)/2.)))/(theta_clonal*((-3*dsample)/2. + (1 - (3*dsample)/2.)*(fbar*phi_clonal + (3*theta_clonal)/2.)))')
    params1.add('c', expr='fbar*phi_clonal/(1+4/3*theta_clonal+fbar*phi_clonal)')
    params1.add('dpool', expr='theta/(1+4/3*theta)')
    params1.add('dclonal', expr='theta_clonal/(1+4/3*theta_clonal)')
    minner1 = Minimizer(fcn2min, params1, fcn_args=(xvalues, yvalues))
    try:
        fitres1 = minner1.minimize()
    except:
        fitres1 = None
    return fitres1

def fit_one(fitdata):
    """Fit one data set"""
    xvalues = fitdata.xvalues
    yvalues = fitdata.yvalues
    dsample = fitdata.d_sample
    fitres = fit_model(xvalues, yvalues, dsample)
    return FitRes(fitdata.group, fitres, dsample)

def fit_p2(fitdatas):
    """Fit p2"""
    all_results = []
    for fitdata in tqdm(fitdatas.getall()):
        fitres = fit_one(fitdata)
        if fitres is not None:
            all_results.append(fitres)
    return all_results


