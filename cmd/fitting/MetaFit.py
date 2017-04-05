#!/usr/bin/env python3
"""Infer recombination rates by fitting mutation correlations"""
from __future__ import print_function
from argparse import ArgumentParser
import numpy as np
from scipy import stats
from lmfit import Parameters, Minimizer
import matplotlib.pyplot as plt

class CorrRes(object):
    """One correlation result"""
    def __init__(self, terms):
        lag = float(terms[0])
        value = float(terms[1])
        variance = float(terms[2])
        num = float(terms[3])
        corrtype = terms[4]
        self.lag = lag
        self.value = value
        self.variance = variance
        self.num = num
        self.corrtype = corrtype

class FitData(object):
    """Fitting data"""
    def __init__(self, xvalues, yvalues, sample_diversity):
        self.xvalues = xvalues
        self.yvalues = yvalues
        self.sample_diversity = sample_diversity

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
    def __str__(self):
        return '{"sample_d": %g, "theta": %g, "phi": %g, "fbar": %g}' % (self.sample_d, self.theta, self.phi, self.fbar)


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

def prepare_fitting_data(corr_results, xmin, xmax):
    """Prepare fitting xvalues and yvalues"""
    p4xs = []
    p4ys = []
    p2xs = []
    p2ys = []
    diver = 0
    for row in corr_results:
        if row.corrtype == "Ks":
            diver = row.value
        elif row.corrtype == "P2" and row.lag > 0 and row.lag >= xmin and row.lag <= xmax:
            p2xs.append(row.lag)
            p2ys.append(row.value)
        elif row.corrtype == "P4" and row.lag > 0 and row.lag >= xmin and row.lag <= xmax:
            p4xs.append(row.lag)
            p4ys.append(row.value)
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(p4ys[:10], p2ys[:10])
    xvalues = []
    yvalues = []
    for i, x in enumerate(p4xs):
        xvalues.append(x)
        y = p4ys[i]
        yvalues.append(y * slope + intercept)
    xvalues = np.array(xvalues)
    yvalues = np.array(yvalues)
    return FitData(xvalues, yvalues, diver)

def calc_p2(diversity, theta, rrate):
    """bulk p2 expression"""
    p2_value = diversity * (1.0 / (2.0 * theta * 4.0 / 3.0 + 4.0 / 3.0 * rrate + 1.0) + 1.0)
    return p2_value

def fcnmin(params, xvalues, yvalues):
    """Model function"""
    theta = params['theta']
    phi = params['phi']
    fbar = params['fbar']
    sample_diversity = params['ds']
    diversity = theta / (1.0 + 4.0 / 3.0 * theta)
    rrate = phi * xvalues
    rcover = sample_diversity / diversity
    if rcover > 1.0:
        rcover = 1.0
    factor = (2.0 * rcover - 1.0 + (1.0 - rcover) ** (1.0 + (1.0 - np.exp(-xvalues/fbar)))) / rcover
    predicts = factor * calc_p2(diversity, theta, rrate)
    return predicts - yvalues

def get_best_model_index(fitresults):
    """Get the index of the best model"""
    aics = []
    for fitres in fitresults:
        aics.append(fitres.aic)
    return np.argmin(aics)

def fit_model1(xvalues, yvalues, sample_diversity):
    """Do fitting using the Model 1"""
    phi_start_values = [0.1]
    fitresults = []
    for phi_start in phi_start_values:
        params1 = Parameters()
        sample_theta = sample_diversity / (1.0 - 4.0/3.0 * sample_diversity)
        params1.add('theta', value=0.1, min=sample_theta)
        params1.add('phi', value=phi_start, min=0, max=1)
        params1.add('fbar', value=100, min=1, max=10000000)
        params1.add('ds', value=sample_diversity, vary=False)
        minner1 = Minimizer(fcnmin, params1, fcn_args=(xvalues, yvalues))
        fitres1 = minner1.minimize()
        fitresults.append(fitres1)
    best_model_index = get_best_model_index(fitresults)
    return fitresults[best_model_index]

def fit_one(fitdata):
    """Fit one data set"""
    fitres = fit_model1(fitdata.xvalues, fitdata.yvalues, fitdata.sample_diversity)
    return (fitdata.group, fitres)

def fit_all(fitdatas):
    """Fit all bootstrap data sets"""
    allresults = Parallel(n_jobs=cpu_count())(delayed(fit_one)(data) for data in fitdatas.getall())
    return allresults

def plot_fit(fitdata, plot_file, json_file):
    """Fit all row data and do ploting"""
    xvalues = fitdata.xvalues
    yvalues = fitdata.yvalues
    sample_diversity = fitdata.sample_diversity
    fitres = fit_model1(xvalues, yvalues, sample_diversity)

    fig = plt.figure(tight_layout=True)
    fig.set_figheight(2)
    fig.set_figwidth(3)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(xvalues, yvalues, 'o',
             markersize=4,
             markeredgewidth=0.75,
             markeredgecolor='k',
             markerfacecolor='None')
    predictions = yvalues + fitres.residual
    ax1.plot(xvalues, predictions, 'k')
    ax1.set_xlabel(r'distance $l$ (bp)')
    ax1.set_ylabel(r'$\tilde P^{(2)}_{s,2}$')
    ax1.locator_params(axis='x', nbins=5)
    ax1.locator_params(axis='y', nbins=5)
    fig.savefig(plot_file)
    with open(json_file, 'w') as writer:
        writer.write(str(FitRes("all", fitres, sample_diversity)))

def fitp2(corr_file, prefix, xmin, xmax):
    """Fit p2"""
    corr_results = read_corr(corr_file)
    fitdata = prepare_fitting_data(corr_results, xmin, xmax)

    plot_file = prefix + "_best_fit.svg"
    json_file = prefix + "_best_fit.json"
    plot_fit(fitdata, plot_file, json_file)

def main():
    """Run fitting using lmfit"""
    parser = ArgumentParser(description="Infer recombination rates\
                                         by fitting mutation correlations.")
    parser.add_argument("corr_file", type=str)
    parser.add_argument("output_prefix", type=str)
    parser.add_argument('--xmin', nargs='?', const=3, type=int, default=3)
    parser.add_argument('--xmax', nargs='?', const=150, type=int, default=150)
    opts = parser.parse_args()
    datafile = opts.corr_file
    prefix = opts.output_prefix
    xmin = opts.xmin
    xmax = opts.xmax

    fitp2(datafile, prefix, xmin, xmax)

if __name__ == "__main__":
    main()
