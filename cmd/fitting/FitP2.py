#!/usr/bin/env python3
"""Infer recombination rates by fitting mutation correlations"""
from __future__ import print_function
from argparse import ArgumentParser
import numpy as np
from lmfit import Parameters, Minimizer
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import cpu_count

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
    def __init__(self, group, xvalues, yvalues, sample_diversity):
        self.group = group
        self.xvalues = xvalues
        self.yvalues = yvalues
        self.sample_diversity = sample_diversity

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
            xvalues, yvalues, sample_diver = prepare_fitting_data(items, xmin, xmax)
            fitdata_map[group] = FitData(group, xvalues, yvalues, sample_diver)
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

def calc_p2(diversity, theta, rrate):
    """bulk p2 expression"""
    p2_value = diversity * (1.0 / (2.0 * theta * 4.0 / 3.0 + 4.0 / 3.0 * rrate + 1.0) + 1.0)
    return p2_value

def constant_one_site(xvalue, fbar):
    """one site function of constant size"""
    return xvalue / fbar

def expon_one_site(xvalue, fbar):
    """one site function of exponential decay size"""
    return 1.0 - np.exp(-xvalue/fbar)

def fcnmin(params, xvalues, yvalues):
    """Model function"""
    global ONESITEFUNC
    theta = params['theta']
    phi = params['phi']
    fbar = params['fbar']
    sample_diversity = params['ds']
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

def fit_model1(xvalues, yvalues, sample_diversity):
    """Do fitting using the Model 1"""
    phi_start_values = [0.001, 0.01, 0.1, 1]
    fitresults = []
    for phi_start in phi_start_values:
        params1 = Parameters()
        sample_theta = sample_diversity / (1.0 - 4.0/3.0 * sample_diversity)
        params1.add('theta', value=0.1, min=sample_theta)
        params1.add('phi', value=phi_start, min=0)
        params1.add('fbar', value=1000, min=3, max=10000000)
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

def plot_fit(fitdata, plot_file):
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

def plot_params(fitresults, param_names, plot_file):
    """Plot histograms"""
    num_col = 3
    num_row = len(param_names) // num_col
    if len(param_names) % num_col > 0:
        num_row = num_row + 1
    fig = plt.figure(tight_layout=True)
    fig.set_figheight(num_row * 6)
    fig.set_figwidth(num_col * 3)
    label_names = {"theta": r"$\theta$",
                   "phi": r"$\phi$",
                   "fbar":r"$\bar f$",
                   "sample_d": r'$d_s$',
                   "ratio": r'$\gamma/\mu$',
                   "rho": r'$\rho$',
                   "sample_theta": r'$\theta_s$',
                   "sample_rho": r'$\rho_s$'}
    for (i, name) in enumerate(param_names):
        values = []
        raw_value = None
        for fitres in fitresults:
            if hasattr(fitres, name):
                value = getattr(fitres, name)
                if "boot" not in fitres.group:
                    raw_value = value
                else:
                    values.append(value)
        if len(values) > 0:
            ax1 = fig.add_subplot(num_col * num_row, num_col, i + 1)
            ax1.hist(values, histtype='stepfilled', bins=30, color="black", alpha=0.5)
            label = label_names.get(name, name)
            ax1.set_xlabel(label)
            ax1.axvline(x=raw_value)
            ax1.locator_params(axis='x', nbins=4)
            ax1.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))

    fig.savefig(plot_file)

def getKey(item):
    """return the first item"""
    return item[0]

def fitp2(corr_file, prefix, xmin, xmax):
    """Fit p2"""
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, xmin, xmax)

    best_fit_plot_file = prefix + "_best_fit.svg"
    plot_fit(fitdatas.get("all"), best_fit_plot_file)

    all_results = sorted(fit_all(fitdatas), key=getKey)
    model_params = ["group", "sample_d", "theta",
                    "phi", "fbar", "ratio", "rho",
                    "sample_theta", "sample_phi", "sample_rho"]
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
    plot_params(fit_results, model_params[1:], out_prefix+".svg")

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
