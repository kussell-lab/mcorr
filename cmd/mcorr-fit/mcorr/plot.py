import numpy
import matplotlib.pyplot as plt
from matplotlib import gridspec

plt.rcParams['mathtext.fontset'] = 'cm'

def plot_fit(fitdata, fitres, plot_file):
    """Fit all row data and do ploting"""
    xvalues = fitdata.xvalues
    yvalues = fitdata.yvalues
    fig = plt.figure(tight_layout=True)

    fig.set_figheight(4)
    fig.set_figwidth(6)
    gs = gridspec.GridSpec(2, 2, 
        height_ratios=[3, 1], width_ratios=[2, 1], hspace=0)
    ax1 = plt.subplot(gs[0, 0])
    ax1.scatter(xvalues, yvalues, s=20, facecolors='none', edgecolors='k')
    predictions = yvalues + fitres.residual
    ax1.plot(xvalues, predictions, 'k')
    ax1.set_ylabel(r'$P$')
    ax1.set_ylim([numpy.min(yvalues)*0.9, numpy.max(yvalues)*1.1])
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

def plot_params(fit_results, param_names, plot_file):
    """plot histogram of parameters"""
    # determine how many columns and rows of sub-plots.
    num_col = 3
    num_row = len(param_names) // num_col
    if len(param_names) % num_col != 0:
        num_row = num_row + 1
    # each sub plot is 6x3
    figheight = num_row * 8
    figwidth = num_col * 3
    # init figure.
    fig = plt.figure(tight_layout=True)
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    # create greek name of parameters.
    label_names = {
            "d_sample": r"$d_{sample}$",
            "fbar": r"$\bar{f}$",
            "theta_pool": r"$\theta_{pool}$",
            "phi_pool": r"$\phi_{pool}$",
            "c": r"$c$",
            "ratio": r"$\gamma/\mu$",
            }

    for i, name in enumerate(param_names):
        boot_values = []
        raw_value = None
        for fitres in fit_results:
            if hasattr(fitres, name):
                value = getattr(fitres, name)
                if fitres.group == "all":
                    raw_value = value
                else:
                    boot_values.append(value)
        if len(boot_values) > 0:
            ax1 = fig.add_subplot(num_col * num_row, num_col, i+1)
            ax1.hist(boot_values, histtype='bar',
                    bins='auto', color="green")
            label = label_names.get(name, name)
            ax1.set_xlabel(label)
            if raw_value is not None:
                ax1.axvline(x=raw_value, color="red")
            ax1.locator_params(axis='x', nbins=6)
            ax1.tick_params(axis='x', rotation=30)
            ax1.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
    fig.savefig(plot_file)

