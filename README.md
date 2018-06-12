# mcorr
Using _Correlation Profile_ of mutations to infer the recombination rate from large-scale sequencing data in bacteria.

## Requirments
* Install `git` from [https://git-scm.com](https://git-scm.com/);
* Install `go` from [https://golang.org/doc/install](https://golang.org/doc/install);
* Install `python3` from [https://www.python.org/](https://www.python.org/) (we found running issues using the default Python in MacOS);
* Install `pip3` from [https://pip.pypa.io/en/stable/installing/](https://pip.pypa.io/en/stable/installing/).

## Installation
1. Install `mcorr-xmfa`, `mcorr-bam`, and `mcorr-fit` from your terminal:
```sh
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-xmfa
go get -u github.com/kussell-lab/mcorr/cmd/mcorr-bam
cd $HOME/go/src/github.com/kussell-lab/mcorr/cmd/mcorr-fit
python3 setup.py install
```
or to install `mcorr-fit` in local directory (~/.local/bin in Linux):
```sh
python3 setup.py install --user
```
2. Add `$HOME/go/bin` and `$HOME/.local/bin` to your `$PATH` environment. In MacOS and Linux, you can do it in your terminal:
```sh
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin
```

We have tested installation in Windows 10, Ubuntu 17.10, and MacOS High Sierra, using Python 3 and Go v1.9.2.

Typical installation time on an iMac is 10 minutes.

## Basic Usage
The inference of recombination parameters requires two steps:

1. Calculate _Correlation Profile_

    For whole-genome alignments (multiple gene alignments), use `mcorr-xmfa`:

    ```sh
    mcorr-xmfa <input XMFA file> <output prefix>
    ```
    The XMFA files should contain only *coding* sequences. The description of XMFA file can be found in [http://darlinglab.org/mauve/user-guide/files.html](http://darlinglab.org/mauve/user-guide/files.html). We provide two useful pipelines to generate whole-genome alignments:
    * from multiple assemblies: [https://github.com/kussell-lab/AssemblyAlignmentGenerator](https://github.com/kussell-lab/AssemblyAlignmentGenerator);
    * from raw reads: [https://github.com/kussell-lab/ReferenceAlignmentGenerator](https://github.com/kussell-lab/ReferenceAlignmentGenerator)

    For read alignments, use `mcorr-bam`:
    ```sh
    mcorr-bam <GFF3 file> <sorted BAM file> <output prefix>
    ```
    The GFF3 file is used for extracting the coding regions of the sorted BAM file.

    Both programs will produce two files:
    * a .csv file stores the calculated Correlation Profile, which will be used for fitting in the next step;
    * a .json file stores the (intermediate) Correlation Profile for each gene.

2. Fit the Correlation Profile using `mcorr-fit`:

    ```sh
    mcorr-fit <.csv file> <output_prefix>
    ```

    It will produce two files:

    * `<output_prefix>_best_fit.svg` shows the plots of the Correlation Profile, fitting, and residuals;
    * `<output_prefix>_fit_reports.txt` shows the summary of the fitted parameters;
    * `<output_prefix>_fit_results.csv` shows the table of fitted parameters;
    * `<output_prefix>_parameter_histograms.svg` shows the distributions of the fitted parameters.

## Examples
1. [Inferring recombination rates of _Helicobacter pylori_ from whole genome sequences of a set of global strains](https://github.com/kussell-lab/Helicobacter_pylori_global_population);
2. [Inferring recombination rates of _Helicobacter pylori_ from reads sequenced from a transformation experiment](https://github.com/kussell-lab/Helicobacter_pylori_transformation_experiments).
