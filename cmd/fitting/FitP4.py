#!/usr/bin/env python3
"""Infer recombination rates by fitting mutation correlations"""
from __future__ import print_function
from argparse import ArgumentParser
from multiprocessing import cpu_count
import numpy as np
from lmfit import Parameters, Minimizer
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import statsmodels.api as sm

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
    def __init__(self, group, xvalues, yvalues, sample_diversity, qfactor):
        self.group = group
        self.xvalues = xvalues
        self.yvalues = yvalues
        self.sample_diversity = sample_diversity
        self.qfactor = qfactor

class FitDatas(object):
    """Fitting data"""
    def __init__(self, corr_results, xmin, xmax):
        corr_map = {}
        for row in corr_results:
            rows = corr_map.get(row.group, [])
            rows.append(row)
            corr_map[row.group] = rows
        fitdata_map = {}
        for group, items in corr_map.items():
            xvalues, yvalues, sample_diver, qfactor = prepare_fitting_data(items, xmin, xmax)
            fitdata_map[group] = FitData(group, xvalues, yvalues, sample_diver, qfactor)
        self.fitdata_dict = fitdata_map
    def get(self, group):
        """return fit data"""
        fitdata = self.fitdata_dict.get(group, None)
        return fitdata
    def getall(self):
        """return all"""
        groups = sorted(self.fitdata_dict.keys())
        return [self.fitdata_dict[group] for group in groups]

class FitRes(object):
    """Fitting results"""
    def __init__(self, group, fit_res, sample_d):
        self.group = group
        self.sample_d = sample_d
        params = fit_res.params.valuesdict()
        if "theta" in params:
            self.theta = params['theta']
        if 'phi' in params:
            self.phi = params['phi']
        if 'fbar' in params:
            self.fbar = params['fbar']
        if 'phi' in params:
            self.ratio = self.phi / self.theta
            if 'fbar' in params:
                self.rho = self.phi * self.fbar
                self.sample_theta = self.sample_d / self.rho
        if 'theta' in params:
            self.sample_rho = self.sample_d / self.theta
        if 'qfactor' in params:
            self.qfactor = params['qfactor']
        sample_rho = getattr(self, 'sample_rho', None)
        fbar = getattr(self, 'fbar', None)
        if sample_rho is not None and fbar is not None:
            self.sample_phi = sample_rho / fbar

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
    p2xvalues = []
    p2yvalues = []
    p4xvalues = []
    p4yvalues = []
    diver = 0
    p4diver = 0
    for row in fitdata:
        if row.corrtype == 'P2' and row.lag >= xmin and row.lag <= xmax:
            p2xvalues.append(row.lag)
            p2yvalues.append(row.value)
        elif row.corrtype == 'Ks':
            diver = row.value
        elif row.corrtype == 'P4' and row.lag == 0:
            p4diver = row.value
        elif row.corrtype == 'P4' and row.lag >= xmin and row.lag <= xmax:
            p4xvalues.append(row.lag)
            p4yvalues.append(row.value)
    linearfitx = p2yvalues[:9]
    linearfity = p4yvalues[:9]
    model = sm.OLS(linearfity, linearfitx)
    results = model.fit()
    slope = results.params[0]
    xvalues = p4xvalues
    yvalues = [y / slope for y in p4yvalues]
    return (xvalues, yvalues, diver, slope)

def calc_p2(diversity, theta, rrate):
    """bulk p2 expression"""
    p2_value = diversity * (1.0 / (2.0 * theta * 4.0 / 3.0 + 4.0 / 3.0 * rrate + 1.0) + 1.0)
    return p2_value

def constant_one_site(xvalue, fbar):
    """one site function of constant size"""
    return xvalue / fbar

def expon_one_site(xvalue, fbar):
    """one site function of exponential decay size"""
    return - np.exp(-np.array(xvalue/fbar)) + 1.0

def fcnmin(params, xvalues, yvalues):
    """Model function"""
    global ONESITEFUNC
    theta = params['theta']
    phi = params['phi']
    fbar = params['fbar']
    sample_diversity = params['ds']
    qfactor = params['qfactor']
    diversity = theta / (1.0 + 4.0 / 3.0 * theta)
    rrate = phi * fbar * ONESITEFUNC(xvalues, fbar)
    rcover = sample_diversity / diversity
    if rcover > 1.0:
        rcover = 1.0
    factor = (2.0 * rcover - 1.0 + (1.0 - rcover) ** (1.0 + ONESITEFUNC(xvalues, fbar))) / rcover
    predicts = factor * calc_p2(diversity, theta, rrate)
    return predicts - yvalues

def get_best_model_index(fitresults):
    """Get the index of the best model"""
    aics = []
    for fitres in fitresults:
        aics.append(fitres.aic)
    return np.argmin(aics)

def fit_model1(xvalues, yvalues, sample_diversity, qfactor):
    """Do fitting using the Model 1"""
    phi_start_values = [0.1]
    fitresults = []
    for phi_start in phi_start_values:
        params1 = Parameters()
        params1.add('theta', value=0.1, min=0, max=1)
        params1.add('phi', value=phi_start, min=0, max=1)
        params1.add('fbar', value=1000, min=3, max=10000000)
        params1.add('ds', value=sample_diversity, vary=False)
        params1.add('qfactor', value=qfactor, vary=False)
        minner1 = Minimizer(fcnmin, params1, fcn_args=(xvalues, yvalues))
        fitres1 = minner1.minimize()
        fitresults.append(fitres1)
    best_model_index = get_best_model_index(fitresults)
    return fitresults[best_model_index]

def fit_one(fitdata):
    """Fit one data set"""
    fitres = fit_model1(fitdata.xvalues, fitdata.yvalues, fitdata.sample_diversity, fitdata.qfactor)
    return (fitdata.group, fitres)

def fit_all(fitdatas):
    """Fit all bootstrap data sets"""
    allresults = Parallel(n_jobs=cpu_count())(delayed(fit_one)(data) for data in fitdatas.getall())
    return allresults

def plot_fit(corr_results, xmin, xmax, plot_file):
    """Fit all row data and do ploting"""
    all_data = FitDatas(corr_results, 3, 300).get('all')
    xvalues = all_data.xvalues
    yvalues = all_data.yvalues
    sample_diversity = all_data.sample_diversity
    qfactor = all_data.qfactor

    fit_data = FitDatas(corr_results, xmin, xmax).get('all')
    fit_xvalues = np.array(fit_data.xvalues)
    fit_yvalues = np.array(fit_data.yvalues)

    rest_xvalues = []
    rest_yvalues = []
    for i in range(len(xvalues)):
        if xvalues[i] > xmax:
            rest_xvalues.append(xvalues[i])
            rest_yvalues.append(yvalues[i])
    fitres = fit_model1(fit_xvalues, fit_yvalues, sample_diversity, qfactor)
    fig = plt.figure(tight_layout=True)
    fig.set_figheight(3)
    fig.set_figwidth(4)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(fit_xvalues, fit_yvalues, 'o',
             markersize=4,
             markeredgewidth=0.75,
             markeredgecolor='k',
             markerfacecolor='None')
    predictions = fit_yvalues + fcnmin(fitres.params, fit_xvalues, fit_yvalues)
    ax1.plot(fit_xvalues, predictions, 'k--')
    ax1.set_xlabel(r'distance $l$ (bp)')
    ax1.set_ylabel(r'$\tilde P^{(4)}_{s,2}$')
    ax1.locator_params(axis='x', nbins=5)
    ax1.locator_params(axis='y', nbins=10)
    fig.savefig(plot_file)

def getKey(item):
    """return the first item"""
    return item[0]

def fitp2(corr_file, prefix, xmin, xmax):
    """Fit p2"""
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, xmin, xmax)

    if fitdatas.get("all"):
        best_fit_plot_file = prefix + "_best_fit.svg"
        plot_fit(corr_results, xmin, xmax, best_fit_plot_file)

    all_results = sorted(fit_all(fitdatas), key=getKey)
    model_params = ["group", "sample_d", "theta",
                    "phi", "fbar", "ratio", "rho",
                    "sample_theta", "sample_phi", "sample_rho", "qfactor"]
    out_prefix = prefix + "_fit_results"
    out_file = out_prefix + ".csv"
    sep = ","
    fit_results = []
    with open(out_file, 'w') as out:
        out.write(sep.join(model_params)+"\n")
        for (group, minres) in all_results:
            sample_d = fitdatas.get(group).sample_diversity
            fit_res = FitRes(group, minres, sample_d)
            values = fit_res.get_values(model_params)
            out.write(sep.join([str(x) for x in values])+"\n")
            fit_results.append(fit_res)

ONESITEFUNC = constant_one_site

def main():
    """Run fitting using lmfit"""
    parser = ArgumentParser(description="Infer recombination rates\
                                         by fitting mutation correlations.")
    parser.add_argument("corr_file", type=str)
    parser.add_argument("output_prefix", type=str)
    parser.add_argument('--xmin', nargs='?', const=3, type=int, default=3)
    parser.add_argument('--xmax', nargs='?', const=300, type=int, default=300)
    parser.add_argument('--onesite', nargs='?', const="const", type=str, default="const")

    opts = parser.parse_args()
    datafile = opts.corr_file
    prefix = opts.output_prefix
    xmin = opts.xmin
    xmax = opts.xmax
    onesite = opts.onesite
    global ONESITEFUNC
    if onesite == "exp":
        ONESITEFUNC = expon_one_site
    else:
        ONESITEFUNC = constant_one_site

    fitp2(datafile, prefix, xmin, xmax)

if __name__ == "__main__":
    main()
