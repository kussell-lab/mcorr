from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from . import fit_p2, read_corr, FitDatas, \
    write_fitting_results, plot_fit, plot_params, write_fitting_reports, \
    geom_r1, const_r1

def main():
    """Run fitting using lmfit, and generate output files and plots"""
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Infer recombination rates\
                    by fitting correlation profile of mutations.")
    parser.add_argument("corr_file", type = str, help='correlation input file')
    parser.add_argument("output_prefix", type=str, help='output file prefix')
    parser.add_argument('--fit_start', type=int, default=3,
                        help='fitting range starts at')
    parser.add_argument('--fit_end', type=int, default=300,
                        help='fitting range ends at')
    parser.add_argument("--use_geom_frag", action="store_true",
                        help='use geometric distribution for fragment sizes')
    parser.add_argument('--quiet', action="store_true")
    parser.add_argument("--title", type=str, help="plot title", default="")
    opts = parser.parse_args()
    corr_file = opts.corr_file
    prefix = opts.output_prefix
    fit_start = opts.fit_start
    fit_end = opts.fit_end
    quiet = opts.quiet
    use_geom_frag = opts.use_geom_frag
    title = opts.title

    # read correlation results and prepare fitting data
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, fit_start, fit_end)
    # do fitting
    r1_func = const_r1
    if use_geom_frag:
        r1_func = geom_r1
    fit_results = fit_p2(fitdatas, r1_func=r1_func, disable_progress_bar=quiet)
    # parameterms to report
    model_params = ["group", "d_sample", "theta_pool",
                    "phi_pool", "ratio", "fbar", "c", "d_pool",
                    "d_clonal", 'theta_clonal', 'phi_clonal']
    # save fitting results into csv file
    csv_file = prefix + "_fit_results.csv"
    write_fitting_results(fit_results, model_params, csv_file)
    # plot the best fit
    best_fit_file = prefix + "_best_fit.svg"
    fitdata = fitdatas.get("all")
    fitres = None
    for res in fit_results:
        if res.group == "all":
            fitres = res
            break
    if fitres is not None:
        plot_fit(fitdata, fitres, best_fit_file, title=title)
    # write fitting report
    report_file = prefix + "_fit_reports.txt"
    write_fitting_reports(fit_results, model_params[1:7], report_file)
    # plot histogram of fitted parameters
    params_hist_file = prefix + "_parameter_histograms.svg"
    plot_params(fit_results, model_params[1:7], params_hist_file)

if __name__ == "__main__":
    main()
