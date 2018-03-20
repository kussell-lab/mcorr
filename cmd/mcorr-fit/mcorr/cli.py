from argparse import ArgumentParser
from . import fit_p2, read_corr, FitDatas, write_fitting_results, plot_fit, plot_params, write_fitting_reports

def main():
    """Run fitting using lmfit"""
    parser = ArgumentParser(description="Infer recombination rates\
                                         by fitting correlation profile of mutations.")
    parser.add_argument("corr_file", type = str)
    parser.add_argument("output_prefix", type=str)
    parser.add_argument('--xmin', nargs='?', const=3, type=int, default=3)
    parser.add_argument('--xmax', nargs='?', const=300, type=int, default=300)
    opts = parser.parse_args()
    corr_file = opts.corr_file
    prefix = opts.output_prefix
    xmin = opts.xmin
    xmax = opts.xmax

    # read correlation results and prepare fitting data
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, xmin, xmax)
    # do fitting
    fit_results = fit_p2(fitdatas)
    # parameterms to report
    model_params = ["group", "d_sample", "theta_pool",
                    "phi_pool", "ratio", "fbar", "c", "d_pool",
                    "d_clonal", 'theta_clonal', 'phi_clonal']
    # save fitting results into csv file
    csv_file = prefix + "_fit_results.csv"
    write_fitting_results(fit_results, model_params, csv_file)
    # plot histogram of fitted parameters
    params_hist_file = prefix + "_parameter_histograms.svg"
    plot_params(fit_results, model_params[1:7], params_hist_file)
    # plot the best fit
    best_fit_file = prefix + "_best_fit.svg"
    fitdata = fitdatas.get("all")
    fitres = None
    for res in fit_results:
        if res.group == "all":
            fitres = res
            break
    if fitres is not None:
        plot_fit(fitdata, fitres, best_fit_file)
    # write fitting report
    report_file = prefix + "_fit_reports.txt"
    write_fitting_reports(fit_results, model_params[1:7], report_file)

if __name__ == "__main__":
    main()
