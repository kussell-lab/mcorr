from . import FitRes
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

