[![poliastro](https://www.poliastro.space/images/logo_text.png)](https://www.poliastro.space/)

Name:         | poliastro
:------------:|:--------------
**Website**:  | <https://www.poliastro.space/>
**Author**:   | Juan Luis Cano Rodríguez [![orcid](https://img.shields.io/badge/id-0000--0002--2187--161X-a6ce39.svg)](http://orcid.org/0000-0002-2187-161X)
**Version**:  | 0.17.0

[![poliastro](https://img.shields.io/circleci/build/gh/poliastro/poliastro/0.17.x?style=flat-square)](https://circleci.com/gh/poliastro/poliastro/?branch=0.17.x)
[![codecov](https://img.shields.io/codecov/c/github/poliastro/poliastro.svg?style=flat-square)](https://codecov.io/github/poliastro/poliastro?branch=0.17.x)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/poliastro/poliastro/0.17.x.svg)](https://results.pre-commit.ci/latest/github/poliastro/poliastro/0.17.x)

[![docs](https://img.shields.io/badge/docs-0.17.x-brightgreen.svg?style=flat-square)](https://docs.poliastro.space/en/0.17.x/?badge=0.17.x)
[![license](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://github.com/poliastro/poliastro/raw/0.17.x/COPYING)
[![doi](https://zenodo.org/badge/11178845.svg?style=flat-square)](https://zenodo.org/badge/latestdoi/11178845)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat-square)](http://www.astropy.org/)
[![mailing](https://img.shields.io/badge/mailing%20list-groups.io-8cbcd1.svg?style=flat-square)](https://groups.io/g/poliastro-dev)
[![Join the chat at http://chat.poliastro.space/](https://img.shields.io/matrix/poliastro:matrix.org.svg?style=flat-square)](http://chat.poliastro.space/)

[![OpenCollective](https://opencollective.com/poliastro/backers/badge.svg)](#backers)
[![OpenCollective](https://opencollective.com/poliastro/sponsors/badge.svg)](#sponsors)

poliastro is an open source (MIT) pure Python library
for interactive Astrodynamics and Orbital Mechanics,
with a focus on ease of use, speed, and quick visualization.
It provides a simple and intuitive API,
and handles physical quantities with units.

Some features include
orbit propagation, solution of the Lambert\'s problem,
conversion between position and velocity vectors and classical orbital elements
and orbit plotting, among others.
It focuses on interplanetary applications,
but can also be used to analyze artificial satellites in Low-Earth Orbit (LEO).

```python
from poliastro.examples import molniya

molniya.plot()
```

![Molniya orbit](https://github.com/poliastro/poliastro/raw/0.17.x/docs/source/examples/molniya.png)

# Documentation

[![docs](https://img.shields.io/badge/docs-0.17.x-brightgreen.svg?style=flat-square)](https://docs.poliastro.space/en/0.17.x/?badge=0.17.x)

Complete documentation, including a user guide and an API reference, can
be read on the wonderful [Read the Docs](https://readthedocs.org/).

<https://docs.poliastro.space/>

# Examples

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/poliastro/poliastro/0.17.x?labpath=index.ipynb)

In the examples directory you can find several Jupyter notebooks with
specific applications of poliastro. You can launch a cloud Jupyter
server using [binder](https://mybinder.org/) to edit the notebooks
without installing anything. Try it out!

<https://mybinder.org/v2/gh/poliastro/poliastro/0.17.x?labpath=index.ipynb>

# Requirements

poliastro requires the following Python packages:

- NumPy, for basic numerical routines
- Astropy, for physical units and time handling
- numba, for accelerating the code
- jplephem, for the planetary ephemerides using SPICE kernels
- matplotlib, for orbit plotting
- plotly, for 2D and 3D interactive orbit plotting
- SciPy, for root finding and numerical propagation

poliastro is supported on Linux, macOS and Windows on Python 3.8 to 3.10.

[![poliastro](https://img.shields.io/circleci/build/gh/poliastro/poliastro/0.17.x?style=flat-square)](https://circleci.com/gh/poliastro/poliastro/?branch=0.17.x)

# Installation

The easiest and fastest way to get the package up and running is to
install poliastro using [conda](http://conda.io):

```bash
    $ conda install poliastro --channel conda-forge
```
Please check out the [documentation for alternative installation
methods](https://docs.poliastro.space/en/0.17.x/installation.html#alternative-installation-methods).

# Problems and suggestions

If for any reason you get an unexpected error message or an incorrect
result, or you want to let the developers know about your use case,
please open a new issue in the [issue
tracker](https://github.com/poliastro/poliastro/issues) and we will try
to answer promptly.

# Contributing

poliastro is a community project, hence all contributions are more than
welcome! For more information, head to
[CONTRIBUTING.md](https://github.com/poliastro/poliastro/blob/0.17.x/CONTRIBUTING.md).

# Support

[![mailing](https://img.shields.io/badge/mailing%20list-groups.io-8cbcd1.svg?style=flat-square)](https://groups.io/g/poliastro-dev)
[![Join the chat at http://chat.poliastro.space/](https://img.shields.io/matrix/poliastro:matrix.org.svg?style=flat-square)](http://chat.poliastro.space/)

Release announcements and general discussion take place on our [Mailing
List](https://groups.io/g/poliastro-dev) .

For further clarifications and discussions, feel free to join Poliastro
[Chat Room](http://chat.poliastro.space/).

# Citing

If you use poliastro on your project, please [drop me a
line](mailto:hello@juanlu.space).

You can also use the DOI to cite it in your publications. This is the
latest one:

[![doi](https://zenodo.org/badge/11178845.svg?style=flat-square)](https://zenodo.org/badge/latestdoi/11178845)

And this is an example citation format:

    Juan Luis Cano Rodríguez et al.. (2015). poliastro: poliastro 0.4.0. Zenodo. 10.5281/zenodo.17462

# License

[![license](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://github.com/poliastro/poliastro/raw/0.17.x/COPYING)

poliastro is released under the MIT license, hence allowing commercial
use of the library. Please refer to the [COPYING](https://github.com/poliastro/poliastro/blob/0.17.x/COPYING) file.

# Credits

## Contributors

This project exists thanks to all the people who contribute!

![Contributors](https://opencollective.com/poliastro/contributors.svg?width=890&button=false)

## Backers

Thank you to all our backers! [Become a backer](https://opencollective.com/poliastro#backer).

[![Backers](https://opencollective.com/poliastro/backers.svg?width=890)](https://opencollective.com/poliastro#backer)

## Sponsors

Support us by becoming a sponsor. Your logo will show up here with a link to your website.
[Become a sponsor](https://opencollective.com/poliastro#sponsor).

[![Sponsors](https://opencollective.com/poliastro/sponsor/0/avatar.svg)](https://opencollective.com/poliastro/sponsor/0/website)

# FAQ

## What's up with the name?

poliastro comes from Polimi, which is the shortened name of the
Politecnico di Milano, the Italian university where I was studying while
writing this software. It's my tiny tribute to a place I came to love.
*Grazie mille!*

## Is poliastro validated?

Yes! poliastro is a community project that strives to be easy to use,
while at the same time producing correct results
[that are validated](https://github.com/poliastro/validation/)
against other [commonly used Astrodynamics software](https://docs.poliastro.space/en/0.17.x/related.html)
such as GMAT and Orekit.

## Can I suggest new features for poliastro?

Sure, we encourage you to [open an issue](https://github.com/poliastro/validation/issues/new)
so we can discuss future feature additions!

## What\'s the future of the project?

poliastro is actively maintained and receiving an influx of new
contributors thanks to the generous sponsorship of Google and the
European Space Agency. The best way to get an idea of the roadmap is to
see the [Milestones](https://github.com/poliastro/poliastro/milestones)
of the project.
