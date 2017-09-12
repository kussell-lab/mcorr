#!/usr/bin/env python3
"""Infer recombination rates by fitting correlation profile"""
from __future__ import print_function
from argparse import ArgumentParser
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
            self.theta = params['theta']
        if 'phi' in params:
            self.phi = params['phi']
        if 'fbar' in params:
            self.fbar = params['fbar']
        if 'phi' in params:
            self.ratio = self.phi / self.theta
            if 'fbar' in params:
                self.rho = self.phi * self.fbar
        if 'c' in params:
            self.c = params['c']
        else:
            self.c = self.d_sample / self.theta

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

def calc_p2_pool(theta, rrate):
    """bulk p2 pool expression"""
    dpool = calc_dpool(theta)
    p2_value = dpool**2 * (1.0 / (2.0 * theta * 4.0 / 3.0 + 4.0 / 3.0 * rrate + 1.0) + 1.0)
    return p2_value

def constant_one_site(xvalue, fbar):
    """one site function of constant size"""
    return xvalue / fbar

def expon_one_site(xvalue, fbar):
    """one site function of exponential decay size"""
    return 1.0 - np.exp(-xvalue/fbar)

def calc_dpool(theta):
    """calculate d_{pool} from \theta_{pool}"""
    return theta/(1.0+4.0/3.0*theta)

def calc_dclonal(dsample, dpool, c):
    """caluclate d_{clonal}"""
    return (dsample - c*dpool)/(1.0 - c)

def calc_c2(c, r1):
    """calculate c2"""
    return 2.0*c - 1.0 + (1.0 - c)**(1.0+r1)

def calc_c1(c, r1):
    """calculate c1"""
    return 2.0 - 2.0*c - 2*(1.0 - c)**(1.0+r1)

def calc_c0(c, r1):
    """calculate c0"""
    return (1.0 - c)**(1.0+r1)

def fcn2min(params, xvalues, yvalues):
    """function 2 min"""
    global ONESITEFUNC
    theta = params['theta']
    phi = params['phi']
    fbar = params['fbar']
    dsample = params['dsample']
    c_value = params['c']
    dpool = params["dpool"]
    dclonal = params["dclonal"]
    r1 = ONESITEFUNC(xvalues, fbar)
    c2 = calc_c2(c_value, r1)
    c1 = calc_c1(c_value, r1)
    c0 = calc_c0(c_value, r1)
    p2pool = calc_p2_pool(theta, phi*fbar*r1)

    p2sample = c2 * p2pool + c1 * dclonal * dpool + c0 * dclonal * dclonal
    predictions = p2sample / dsample
    return predictions - yvalues

def fit_model(xvalues, yvalues, d_sample, c_start):
    """Do fitting using the Model 1"""
    params1 = Parameters()
    params1.add('dsample', value=d_sample, vary=False)
    params1.add('c', value=c_start, min=0, max=1)
    params1.add('fbar', value=100, min=3, max=10000000)
    params1.add('dclonal', value=d_sample*0.001, min=0, max=d_sample)
    params1.add('dpool', expr="(dsample-(1-c)*dclonal)/c")
    params1.add('theta', expr="dpool/(1-4.0/3.0*dpool)")
    params1.add('theta_clonal', expr="dclonal/(1-4.0/3.0*dclonal)")
    params1.add('phi', expr="-theta*log(1-c)/(theta_clonal*fbar)")
    minner1 = Minimizer(fcn2min, params1, fcn_args=(xvalues, yvalues))
    fitres1 = minner1.minimize()
    return fitres1

def fit_one(fitdata, c_start):
    """Fit one data set"""
    xvalues = fitdata.xvalues
    yvalues = fitdata.yvalues
    dsample = fitdata.d_sample
    fitres = fit_model(xvalues, yvalues, dsample, c_start)
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

def fitp2(corr_file, prefix, xmin, xmax, fit_bootstraps=False, c_start=0.1):
    """Fit p2"""
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, xmin, xmax)

    all_results = []
    if fitdatas.has("all"):
        fitdata = fitdatas.get("all")
        best_fit_plot_file = prefix + "_best_fit.svg"
        fitres = fit_one(fitdata, c_start)
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
            fitres = fit_one(fitdata, c_start)
            all_results.append((fitdata.group, fitres))

    # write fitting results.
    model_params = ["group", "d_sample", "theta",
                    "phi", "fbar", "c", "ratio"]
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

ONESITEFUNC = constant_one_site

def main():
    """Run fitting using lmfit"""
    parser = ArgumentParser(description="Infer recombination rates\
                                         by fitting mutation correlations.")
    parser.add_argument("corr_file", type=str)
    parser.add_argument("output_prefix", type=str)
    parser.add_argument('--xmin', nargs='?', const=3, type=int, default=3)
    parser.add_argument('--xmax', nargs='?', const=150, type=int, default=150)
    parser.add_argument('--onesite', nargs='?', const="const", type=str, default="const")
    parser.add_argument('--fit_bootstraps', nargs='?', const="false", type=bool, default=False)
    parser.add_argument('--c_start', nargs='?', const=0.1, type=float, default=0.1)
    opts = parser.parse_args()
    datafile = opts.corr_file
    prefix = opts.output_prefix
    xmin = opts.xmin
    xmax = opts.xmax
    onesite = opts.onesite
    fit_bootstraps = opts.fit_bootstraps
    c_start = opts.c_start
    global ONESITEFUNC
    if onesite == "exp":
        ONESITEFUNC = expon_one_site
    else:
        ONESITEFUNC = constant_one_site

    fitp2(datafile, prefix, xmin, xmax, fit_bootstraps, c_start)

if __name__ == "__main__":
    main()
