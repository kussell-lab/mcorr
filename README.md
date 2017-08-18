# mcorr
Infer recombination rates from bacterial sequence data using correlated mutations.

## Requirments
* [Git](https://git-scm.com/);
* [Go](https://golang.org/);
* Python with numpy, lmfit, tqdm and matplotlib packages.

## Installation
```sh
go get -u github.com/mingzhi/mcorr/cmd/mcorr-xmfa
go get -u github.com/mingzhi/mcorr/cmd/mcorr-bam
```

## Usage
The inference requires two steps:

1. calculate mutation correlation using `mcorr-xmfa` from gene alignments or `mcorr-bam` from read alignments:

    ```sh
    mcorr-xmfa <input XMFA file)> <output CSV file>
    ```
    or
    ```sh
    mcorr-bam <GFF3 file> <input BAM file> <output (CSV file)>
    ```

    **Note**: the XMFA files should contain only *coding* sequences.
    The description of XMFA file can be found [here](http://darlinglab.org/mauve/user-guide/files.html).

2. fit the correlation results using `FitP2.py`:

    ```sh
    python $GOPATH/src/github.com/mingzhi/mcorr/cmd/fitting/FitP2.py <input (mcorr output file)> <output prefix>
    ```

    The resulted files:

    * `<output_prefix>_best_fit.svg` -- the correlation profile, fitting, and residual plots;
    * `<output_prefix>_fit_results.csv` -- table of the best-fitted parameters.

Example data can be found [here](https://github.com/mingzhi/mcorr_examples).
