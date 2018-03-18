#!/usr/bin/env python3
"""Infer recombination rates by fitting correlation profile"""
from __future__ import print_function
import numpy as np
from lmfit import Parameters, Minimizer
import matplotlib.pyplot as plt
from matplotlib import gridspec
from tqdm import tqdm

plt.rcParams['mathtext.fontset'] = 'cm'

class CorrRes(object):
    """One correlation result"""
    def __init__(self, terms):
        lag = float(terms[0])
        value = float(terms[1])
        variance = float(terms[2])
        num = float(terms[3])
        corrtype = terms[4]
        group = terms[5]
        self.lag = lag
        self.value = value
        self.variance = variance
        self.num = num
        self.corrtype = corrtype
        self.group = group

class FitData(object):
    """Fitting data"""
    def __init__(self, group, xvalues, yvalues, d_sample):
        self.group = group
        self.xvalues = xvalues
        self.yvalues = yvalues
        self.d_sample = d_sample

class FitDatas(object):
    """Fitting data"""
    def __init__(self, corr_results, xmin, xmax):
        corr_map = {}
        groups = []
        for row in corr_results:
            rows = corr_map.get(row.group, [])
            if len(rows) == 0:
                groups.append(row.group)
            rows.append(row)
            corr_map[row.group] = rows
        fitdata_map = {}
        for group, items in corr_map.items():
            xvalues, yvalues, d_sample = prepare_fitting_data(items, xmin, xmax)
            fitdata_map[group] = FitData(group, xvalues, yvalues, d_sample)
        self.fitdata_dict = fitdata_map
        self.groups = groups
    def has(self, group):
        """return True if the group is in the data"""
        return group in self.fitdata_dict

    def get(self, group):
        """return fit data"""
        fitdata = self.fitdata_dict.get(group, None)
        return fitdata
    def getall(self):
        """return all"""
        return [self.fitdata_dict[group] for group in self.groups]

class FitRes(object):
    """Fitting results"""
    def __init__(self, group, fit_res, d_sample):
        self.group = group
        self.d_sample = d_sample
        params = fit_res.params.valuesdict()
        if "theta" in params:
            self.theta_pool = params['theta']
        if 'phi' in params:
            self.phi_pool = params['phi']
        if 'fbar' in params:
            self.fbar = params['fbar']
        if 'phi' in params:
            self.ratio = self.phi_pool / self.theta_pool
            if 'fbar' in params:
                self.rho = self.phi_pool * self.fbar
        if 'c' in params:
            self.c = params['c']
        if 'dclonal' in params:
            self.d_clonal = params['dclonal']
        if 'dpool' in params:
            self.d_pool = params['dpool']
        if 'phi_clonal' in params:
            self.phi_clonal = params['phi_clonal']
        if 'theta_clonal' in params:
            self.theta_clonal = params['theta_clonal']

    def get_values(self, attributes):
        """Get attribute values"""
        values = []
        for name in attributes:
            if hasattr(self, name):
                values.append(getattr(self, name))
            else:
                values.append("NA")
        return values

def read_corr(csv_file):
    """Read corr results in a csv file"""
    results = []
    with open(csv_file, 'r') as infile:
        for line in infile:
            terms = line.rstrip().split(",")
            if terms[0] == 'l':
                continue
            results.append(CorrRes(terms))
    return results

def prepare_fitting_data(fitdata, xmin, xmax):
    """Prepare fitting xvalues and yvalues"""
    xvalues = []
    yvalues = []
    diver = 0
    for row in fitdata:
        if row.corrtype == 'P2' and row.lag >= xmin and row.lag <= xmax:
            xvalues.append(row.lag)
            yvalues.append(row.value)
        elif row.corrtype == 'Ks':
            diver = row.value
    xvalues = np.array(xvalues)
    yvalues = np.array(yvalues)
    return (xvalues, yvalues, diver)

def Power(a, b):
    """compute power"""
    return a**b

def const_r1(x, fBar, phiC):
    """calculate r1 assuming constant fragment size"""
    return np.where(x < fBar, phiC*x, phiC*fBar)

def exp_r1(x, fBar, phiC):
    """calculate r1 assuming exponetional decay of fragment size"""
    return phiC*fBar*(1.0 - np.exp(-x/fBar))

def geom_r1(x, fBar, phiC):
    """calculate r1 assuming geom distribution"""
    prob = 1.0/fBar
    return phiC*fBar*(1.0 - np.power(1-prob, x))

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
    return fitres

def plot_fit(fitdata, fitres, plot_file):
    """Fit all row data and do ploting"""
    xvalues = fitdata.xvalues
    yvalues = fitdata.yvalues
    fig = plt.figure(tight_layout=True)

    fig.set_figheight(4)
    fig.set_figwidth(6)
    gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1], width_ratios=[2, 1], hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1.scatter(xvalues, yvalues, s=20, facecolors='none', edgecolors='k')
    predictions = yvalues + fitres.residual
    ax1.plot(xvalues, predictions, 'k')
    ax1.set_ylabel(r'$P$')
    ax1.set_ylim([np.min(yvalues)*0.9, np.max(yvalues)*1.1])
    ax1.locator_params(axis='x', nbins=5)
    ax1.locator_params(axis='y', nbins=5)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2 = plt.subplot(gs[1, 0])
    markerline, _, _ = ax2.stem(xvalues,
                                fitres.residual,
                                linefmt='k-',
                                basefmt='r-',
                                markerfmt='ko')
    ax2.set_xlabel(r'$l$')
    ax2.set_ylabel("Residual")
    ax2.locator_params(axis='x', nbins=5)
    plt.setp(markerline, "markersize", 4)
    fig.tight_layout()

    ax3 = plt.subplot(gs[1, 1])
    ax3.hist(fitres.residual, bins="auto", facecolor='green', alpha=0.5)
    ax3.set_xlabel("Residual")
    plt.setp(ax3.get_xticklabels(), rotation=10, horizontalalignment='right')
    ax3.axes.get_yaxis().set_ticks([])
    fig.savefig(plot_file)

def fitp2(corr_file, prefix, xmin, xmax, fit_bootstraps=False):
    """Fit p2"""
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, xmin, xmax)

    all_results = []
    if fitdatas.has("all"):
        fitdata = fitdatas.get("all")
        best_fit_plot_file = prefix + "_best_fit.svg"
        fitres = fit_one(fitdata)
        if fitres:
            plot_fit(fitdata, fitres, best_fit_plot_file)
            all_results.append((fitdata.group, fitres))

    to_fit_groups = []
    for fitdata in fitdatas.getall():
        tofit = True
        if fitdata.group == "all":
            tofit = False
        if "boot" in fitdata.group:
            tofit = fit_bootstraps
        if tofit:
            to_fit_groups.append(fitdata.group)
    num_groups = len(to_fit_groups)
    if num_groups > 0:
        for group in tqdm(to_fit_groups):
            fitdata = fitdatas.get(group)
            fitres = fit_one(fitdata)
            if fitres:
                all_results.append((fitdata.group, fitres))

    # write fitting results.
    model_params = ["group", "d_sample", "theta_pool",
                    "phi_pool", "ratio", "fbar", "c", "d_pool",
                    "d_clonal", 'theta_clonal', 'phi_clonal']
    out_prefix = prefix + "_fit_results"
    out_file = out_prefix + ".csv"
    sep = ","
    with open(out_file, 'w') as out:
        out.write(sep.join(model_params)+"\n")
        for (group, minres) in all_results:
            d_sample = fitdatas.get(group).d_sample
            fit_res = FitRes(group, minres, d_sample)
            values = fit_res.get_values(model_params)
            out.write(sep.join([str(x) for x in values])+"\n")

