"""Fit correlation function (P2)"""
from __future__ import print_function
from multiprocessing import cpu_count
from argparse import ArgumentParser
import numpy as np
from joblib import Parallel, delayed
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 12})
import matplotlib.pyplot as plt
from lmfit import Parameters, Minimizer

class Row(object):
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

def read_csv(filename):
    """Read cov results"""
    results = []
    with open(filename, 'r') as infile:
        for line in infile:
            terms = line.rstrip().split(",")
            if terms[0] == 'l':
                continue
            results.append(Row(terms))
    return results

def fcn1min(params, xvalues, yvalues):
    """Model 1: assume the sample is closely related"""
    theta = params['theta']
    phi = params['phi']
    fbar = params['fbar']
    diversity = theta / (1.0 + 4.0 / 3.0 * theta)
    rrate = phi * xvalues
    predicts = diversity * (1.0 / (2.0 * theta * 4.0 / 3.0 + 4.0 / 3.0 * rrate + 1.0) + 1.0) \
     * (fbar - xvalues) / fbar
    return predicts - yvalues

def fcn2min(params, xvalues, yvalues):
    """Model 2: assume the sample is diverse"""
    theta = params['theta']
    phi = params['phi']
    diversity = theta / (1.0 + 4.0 / 3.0 * theta)
    rrate = phi * xvalues
    predicts = diversity * (1.0 / (2.0 * theta * 4.0 / 3.0 + 4.0 / 3.0 * rrate + 1.0) + 1.0)
    return predicts - yvalues

def fcn3min(params, xvalues, yvalues):
    """Model 3: assume there is no recombination"""
    theta = params['theta']
    fbar = params['fbar']
    diversity = theta / (1.0 + 4.0 / 3.0 * theta)
    predicts = diversity * (1.0 / (2.0 * theta * 4.0 / 3.0 + 1.0) + 1.0) * (fbar - xvalues) / fbar
    return predicts - yvalues

def fit_model1(xvalues, yvalues):
    """Do fitting using the Model 1"""
    phi_start_values = [0.001, 0.01, 0.1, 1]
    fitresults = []
    for phi_start in phi_start_values:
        params1 = Parameters()
        params1.add('theta', value=0.1, min=0)
        params1.add('phi', value=phi_start, min=0)
        params1.add('fbar', value=100, min=1, max=10000000)
        minner1 = Minimizer(fcn1min, params1, fcn_args=(xvalues, yvalues))
        fitres1 = minner1.minimize()
        fitresults.append(fitres1)
    best_model_index = get_best_model_index(fitresults)
    return fitresults[best_model_index]

def fit_model2(xvalues, yvalues):
    """Do fitting using the Model 2"""
    phi_start_values = [0.001, 0.01, 0.1, 1]
    fitresults = []
    for phi_start in phi_start_values:
        params2 = Parameters()
        params2.add('theta', value=0.1, min=0)
        params2.add('phi', value=phi_start, min=0)
        minner2 = Minimizer(fcn2min, params2, fcn_args=(xvalues, yvalues))
        fitres2 = minner2.minimize()
        fitresults.append(fitres2)
    best_model_index = get_best_model_index(fitresults)
    return fitresults[best_model_index]

def fit_model3(xvalues, yvalues):
    """Do fitting using the Model 3"""
    params3 = Parameters()
    params3.add('theta', value=0.1, min=0)
    params3.add('fbar', value=100, min=1)
    minner3 = Minimizer(fcn3min, params3, fcn_args=(xvalues, yvalues))
    fitres3 = minner3.minimize()
    return fitres3

def fit(xvalues, yvalues):
    """Fit three models for one fitting data"""
    fitres1 = fit_model1(xvalues, yvalues)
    fitres2 = fit_model2(xvalues, yvalues)
    fitres3 = fit_model3(xvalues, yvalues)
    return (fitres1, fitres2, fitres3)

def fit_one(data):
    """Fit only one"""
    group, fitdata = data
    xvalues, yvalues, _ = prepare_fitting_data(fitdata)
    fitres1, fitres2, fitres3 = fit(xvalues, yvalues)
    results = [fitres1, fitres2, fitres3]
    return (group, results)

def fit_all(fitdatas):
    """Fit all bootstrap data sets"""
    allresults = Parallel(n_jobs=cpu_count())(delayed(fit_one)(data) for data in fitdatas.items())
    return allresults

def plot_fit(fitdatas, ax1):
    """Fit all row data and do ploting"""
    fitdata = fitdatas['all']
    xvalues, yvalues, yerrors = prepare_fitting_data(fitdata)
    fitres1, fitres2, fitres3 = fit(xvalues, yvalues)
    fitres_dict = {"Model_1": fitres1, "Model_2": fitres2, "Model_3": fitres3}
    color_dict = {"Model_1": 'b', "Model_2": 'r', "Model_3": 'g'}
    ax1.errorbar(xvalues, yvalues, yerr=yerrors, fmt="ko")
    for name, fitres in fitres_dict.items():
        predictions = yvalues + fitres.residual
        ax1.plot(xvalues, predictions, color_dict[name], label=name, lw=2)
    ax1.legend(loc='upper right')
    ax1.set_xlabel(r'distance $l$ (bp)')
    ax1.set_ylabel(r'$\tilde P^{(2)}_{s,2}$')
    ax1.locator_params(axis='x', nbins=4)
    ax1.grid(True)
    ax1.set_title("Best fitting results for different models")

def prepare_fitting_data(fitdata):
    """Prepare fitting xvalues and yvalues"""
    xvalues = []
    yvalues = []
    yerrors = []
    for row in fitdata:
        xvalues.append(row.lag)
        yvalues.append(row.value)
        yerrors.append(np.sqrt(row.variance/row.num))
    xvalues = np.array(xvalues)
    yvalues = np.array(yvalues)
    return (xvalues, yvalues, yerrors)

def get_best_model_index(fitresults):
    """Get the index of the best model"""
    aics = []
    for fitres in fitresults:
        aics.append(fitres.aic)
    return np.argmin(aics)

def plot_aics(allresults, ax1):
    """Plot histogram of AICs and BICs"""
    aics = [[], [], []]
    fit_aics = []
    for group, results in allresults:
        for i in range(3):
            res = results[i]
            if "boot" not in group:
                fit_aics.append(res.aic)
            else:
                aics[i].append(res.aic)
    name = "AIC"
    ax1.hist(aics[0], label="model_1", alpha=0.5, histtype='stepfilled', bins=30, color='b')
    ax1.hist(aics[1], label="model_2", alpha=0.5, histtype='stepfilled', bins=30, color='r')
    ax1.hist(aics[2], label="model_3", alpha=0.5, histtype='stepfilled', bins=30, color='g')
    ax1.axvline(x=fit_aics[0], color='b')
    ax1.axvline(x=fit_aics[1], color='r')
    ax1.axvline(x=fit_aics[2], color='g')
    ax1.legend(loc='upper right')
    ax1.grid(True)
    ax1.set_xlabel(name)
    ax1.set_title("Comparison of AIC for different models\n in bootstrap results")

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
            ax1.grid(True)
            ax1.locator_params(axis='x', nbins=4)
            ax1.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))

    fig.savefig(plot_file)


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
    def get_values(self, attributes):
        """Get attribute values"""
        values = []
        for name in attributes:
            if hasattr(self, name):
                values.append(getattr(self, name))
            else:
                values.append("NA")
        return values

def main():
    """Run fitting using lmfit"""
    parser = ArgumentParser(description="fitting p2")
    parser.add_argument("corr_file", type=str)
    parser.add_argument("out_prefix", type=str)
    opts = parser.parse_args()
    datafile = opts.corr_file
    prefix = opts.out_prefix
    data = read_csv(datafile)
    xmin = 3
    xmax = 150
    fitdatas = {}
    sample_d_dict = {}
    for row in data:
        if row.corrtype == 'P2' and row.lag >= xmin and row.lag <= xmax:
            group = row.group
            if group not in fitdatas:
                fitdatas[group] = []
            fitdatas[group].append(row)
        elif row.corrtype == 'Ks' and row.lag == 0:
            sample_d_dict[row.group] = row.value
    all_results = fit_all(fitdatas)
    _, (ax1, ax2) = plt.subplots(2, 1)
    plot_fit(fitdatas, ax1)
    plot_aics(all_results, ax2)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(6, 8)
    plt.tight_layout()
    fig.savefig(prefix + "_fit.svg")

    model1_params = ["group", "sample_d", "theta",
                     "phi", "fbar", "ratio", "rho",
                     "sample_theta", "sample_rho"]
    model2_params = ["group", "sample_d", "theta",
                     "phi", "ratio", "sample_rho"]
    model3_params = ["group", "sample_d", "theta",
                     "fbar", "sample_rho"]
    param_names = [model1_params, model2_params, model3_params]
    for index in range(3):
        out_prefix = prefix + "_model_%d_fit_results" % (index + 1)
        out_file = out_prefix + ".csv"
        names = param_names[index]
        sep = ","
        fit_results = []
        with open(out_file, 'w') as out:
            out.write(sep.join(names)+"\n")
            for (group, results) in all_results:
                sample_d = sample_d_dict[group]
                fit_res = FitRes(group, results[index], sample_d)
                values = fit_res.get_values(names)
                out.write(sep.join([str(x) for x in values])+"\n")
                fit_results.append(fit_res)
        plot_params(fit_results, names[1:], out_prefix+".svg")

if __name__ == "__main__":
    main()
