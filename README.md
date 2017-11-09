# mcorr
Using correlation profile of mutations to infer bacterial recombination rates from large-scale sequencing datasets.

## Requirments
* Install `git` from [https://git-scm.com](https://git-scm.com/);
* Install `go` from [https://golang.org](https://golang.org/);
* Install `python` from [https://www.python.org/](https://www.python.org/) and `pip` from [https://pypi.python.org/pypi/pip/](https://pypi.python.org/pypi/pip/).
* Use `pip` to install `numpy`, `matplotlib`, `lmfit`, and `tqdm`:

    `pip install --user numpy matplotlib lmfit tqdm`

## Installation
1. Download and install `mcorr-xmfa` and `mcorr-bam` from your terminal:
```sh
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-xmfa
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-bam
```
2. Add `$HOME/go/bin` to your `$PATH` environment. Both programs are installed in `$HOME/go/bin`.

## Usage
The inference requires two steps:

1. calculate correlation profile of mutations using `mcorr-xmfa` for gene alignments:

    ```sh
    mcorr-xmfa <input XMFA file> <output prefix>
    ```
    The XMFA files should contain only *coding* sequences. The description of XMFA file can be found in [http://darlinglab.org/mauve/user-guide/files.html](http://darlinglab.org/mauve/user-guide/files.html).

    Or use `mcorr-bam` for raw read alignments:
    ```sh
    mcorr-bam <GFF3 file> <sorted BAM file> <output prefix>
    ```
    The GFF3 file is used for extracting the coding regions of the sorted BAM file.

    Both programs will produce two files:
    * a <output prefix>.csv file stores the computed correlation profile, which will be used for fitting in the next step;
    * a <output prefix>.json file stores the (intermediate) correlation results for each gene.

2. fit the correlation profile using `FitP.py`, which is located in `$HOME/go/src/github.com/kussell-lab/mcorr/cmd/fitting/`:

    ```sh
    python $HOME/go/src/github.com/kussell-lab/mcorr/cmd/fitting/FitP.py <input (the .csv file)> <output prefix>
    ```

    It will produce two files:

    * `<output_prefix>_best_fit.jpg` -- the plots of correlation profile, fitting, and residuals;
    * `<output_prefix>_fit_results.csv` -- the table of fitted parameters.

Example data can be found [here](https://github.com/kussell-lab/mcorr_examples).
