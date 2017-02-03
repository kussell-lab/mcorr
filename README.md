# mcorr
Infer recombination rates from bacterial sequence data using correlated mutations.

## Requirments
* [Go](https://golang.org/) (be sure to setup GOPATH).
* Python with numpy, lmfit and matplotlib.
* [Git](https://git-scm.com/)

## Installation
```sh
go get -u github.com/mingzhi/mcorr/cmd/mcorr
```

## Usage
The inference of recombination rates requires two steps:

1. calculate mutation correlation from sequence data using `mcorr`;

    ```sh
    mcorr <input ([XMFA] file)> <output (CSV file)>
    ```

    **Note**: the program assumes that the input file contains only *coding* sequences.
    The description of XMFA file can be found in [here](http://darlinglab.org/mauve/user-guide/files.html).

2. fit the output correlation results using `FitP2.py`:

    ```sh
    python $GOPATH/src/github.com/mingzhi/mcorr/cmd/fitting/FitP2.py <input (mcorr output file)> <output prefix>
    ```

    It produces three files:

    * `<output_prefix>_best_fit.svg`
    * `<output_prefix>_fit_results.csv` -- table of the best-fitted parameters
    * `<output_prefix>_fit_results.svg` -- histograms of the best-fitted parameters