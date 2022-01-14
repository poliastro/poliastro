# Installation

## Requirements

poliastro requires a number of Python packages, notably:

- [Astropy](https://www.astropy.org/), for physical units and time handling
- [NumPy](https://numpy.org/), for basic numerical routines
- [jplephem](https://pypi.org/project/jplephem/), for the planetary ephemerides using SPICE kernels
- [matplotlib](https://matplotlib.org/), for static orbit plotting
- [numba](https://numba.pydata.org/) (when using CPython), for accelerating the code
- [Plotly](https://plotly.com/), for interactive orbit plotting
- [SciPy](https://www.scipy.org/), for root finding and numerical propagation

poliastro is supported on Linux, macOS and Windows on Python 3.8 to 3.9.

## Using conda

The easiest and fastest way to get the package up and running is to
install poliastro using [conda](https://conda.io/docs/):

```bash
$ conda install -c conda-forge poliastro
```

or, better yet, using [mamba](https://mamba.readthedocs.io/),
which is a super fast replacement for `conda`:

```bash
$ conda install -c conda-forge mamba
$ mamba install -c conda-forge poliastro
```

```{note}
We encourage users to use conda or mamba
and the [conda-forge](https://conda-forge.org/) packages
for convenience,
especially when developing on Windows.
It is recommended to create a new environment.
```

If the installation fails for any reason, please open an issue in the
[issue tracker](https://github.com/poliastro/poliastro/issues).

## Alternative installation methods

You can also [install poliastro from PyPI](https://pypi.python.org/pypi/poliastro/) using pip:

```bash
$ pip install poliastro
```

Finally, you can also install the latest development version of poliastro
[directly from GitHub](http://github.com/poliastro/poliastro):

```bash
$ pip install https://github.com/poliastro/poliastro/archive/main.zip
```

This is useful if there is some feature that you want to try,
but we did not release it yet as a stable version.
Although you might find some unpolished details,
these development installations should work without problems.
If you find any, please open an issue in the [issue tracker](https://github.com/poliastro/poliastro/issues).

```{warning}
It is recommended that you
**never ever use sudo** with distutils, pip, setuptools and friends in Linux
because you might seriously break your system
\[[1](http://wiki.python.org/moin/CheeseShopTutorial#Distutils_Installation)\]\[[2](http://stackoverflow.com/questions/4314376/how-can-i-install-a-python-egg-file/4314446#comment4690673_4314446)\]\[[3](http://workaround.org/easy-install-debian)\]\[[4](http://matplotlib.1069221.n5.nabble.com/Why-is-pip-not-mentioned-in-the-Installation-Documentation-tp39779p39812.html)\].
Use [virtual environments](https://docs.python.org/3/library/venv.html) instead.
```

## Making poliastro work in your editor

### Jupyter notebook and JupyterLab

To install the extra dependencies needed to make the interactive plots work on Jupyter, do

```bash
$ pip install poliastro[jupyter]
```

With Plotly versions older than 5 on JupyterLab,
you will also need to install Node.js
to enable the browser extensions.
Check out [their troubleshooting guide](https://plotly.com/python/troubleshooting/#jupyterlab-problems)
for further information.

## Problems and suggestions

If for any reason you get an unexpected error message or an incorrect result,
or you want to let the developers know about your use case,
please open a new issue in the [issue tracker](https://github.com/poliastro/poliastro/issues)
and we will try to answer promptly.
