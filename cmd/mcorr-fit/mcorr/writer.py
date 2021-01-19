from . import FitReport
def write_fitting_results(all_results, model_params, out_file):
    """
    write fitting results into a .csv file.
    """
    # write fitting results.
    sep = ","
    with open(out_file, 'w') as out:
        out.write(sep.join(model_params)+"\n")
        for fit_res in all_results:
            values = fit_res.get_values(model_params)
            out.write(sep.join([str(x) for x in values])+"\n")

def write_fitting_reports(all_results, model_params, out_file):
    """
    write fitting reports into a .txt file.
    """
    with open(out_file, 'w') as out:
        for param_name in model_params:
            label_name = param_name
            if param_name == "ratio":
                label_name = "gamma/mu"
            report = FitReport(all_results, param_name, label_name)
            out.write(report.report()+"\n")

