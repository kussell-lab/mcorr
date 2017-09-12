# mcorr
Inferring bacterial recombination rates from large-scale sequencing datasets.

## Requirments
* [Git](https://git-scm.com/);
* [Go](https://golang.org/);
* Python with numpy, lmfit, tqdm and matplotlib packages.

## Installation
```sh
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-xmfa
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-bam
```

## Usage
The inference requires two steps:

1. calculate mutation correlation using `mcorr-xmfa` from gene alignments or `mcorr-bam` from read alignments:

    ```sh
    mcorr-xmfa <input XMFA file)> <output prefix>
    ```
    or
    ```sh
    mcorr-bam <GFF3 file> <input BAM file> <output prefix>
    ```

    **Note**: the XMFA files should contain only *coding* sequences.
    The description of XMFA file can be found [here](http://darlinglab.org/mauve/user-guide/files.html).

2. fit the correlation results using `FitP.py`:

    ```sh
    python $HOME/go/src/github.com/kussell-lab/mcorr/cmd/fitting/FitP.py <input (mcorr output file)> <output prefix>
    ```

    The resulted files:

    * `<output_prefix>_best_fit.svg` -- the correlation profile, fitting, and residual plots;
    * `<output_prefix>_fit_results.csv` -- table of the best-fitted parameters.

Example data can be found [here](https://github.com/kussell-lab/mcorr_examples).
