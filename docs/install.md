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
pip3 install --user mcorr
```
2. Add `$HOME/go/bin` and `$HOME/.local/bin` to your `$PATH` environment:
```sh
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin
```

We have tested installation in Windows 10, Ubuntu 17.10, and MacOS High Sierra, using Python 3 and Go v1.9.2.

Typical installation time on an iMac is 10 minutes.

