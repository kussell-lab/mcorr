# mcorr
Inferring bacterial recombination rates from large-scale sequencing datasets.

## Requirments
* Install `git` from [https://git-scm.com](https://git-scm.com/);
* Install `go` from [https://golang.org](https://golang.org/);
* Install `python` from [https://www.python.org/](https://www.python.org/) and `pip` from [https://pypi.python.org/pypi/pip/](https://pypi.python.org/pypi/pip/).
* Use `pip` to install `numpy`, `matplotlib`, `lmfit`, and `tqdm`:

    `pip install --user numpy matplotlib lmfit tqdm`

## Installation
```sh
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-xmfa
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-bam
```

## Usage
The inference requires two steps:

1. calculate mutation correlation using `mcorr-xmfa` for gene alignments or `mcorr-bam` for raw read alignments:

    ```sh
    mcorr-xmfa <input XMFA file> <output prefix>
    ```
    or
    ```sh
    mcorr-bam <GFF3 file> <input BAM file> <output prefix>
    ```

    **Note**: the XMFA files should contain only *coding* sequences.
    The description of XMFA file can be found in [http://darlinglab.org/mauve/user-guide/files.html](http://darlinglab.org/mauve/user-guide/files.html).
    The GFF3 file is the .gff3 file of the reference genome used for read mapping.

2. fit the correlation results using `FitP.py`:

    ```sh
    python $HOME/go/src/github.com/kussell-lab/mcorr/cmd/fitting/FitP.py <input (mcorr output file)> <output prefix>
    ```

    The resulted files:

    * `<output_prefix>_best_fit.jpg` -- the plots of correlation profile, fitting, and residuals;
    * `<output_prefix>_fit_results.csv` -- the table of fitted parameters.

Example data can be found [here](https://github.com/kussell-lab/mcorr_examples).
