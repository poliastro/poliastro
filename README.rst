.. poliastro

.. image:: http://poliastro.github.io/images/logo_text.png
   :target: http://poliastro.github.io/
   :alt: poliastro logo
   :width: 675px
   :align: center

.. |orcid| image:: https://img.shields.io/badge/id-0000--0002--2187--161X-a6ce39.svg
   :target: http://orcid.org/0000-0002-2187-161X

:Name: poliastro
:Website: https://poliastro.github.io/
:Author: Juan Luis Cano Rodríguez |orcid|
:Version: 0.13.dev0

.. |circleci| image:: https://img.shields.io/circleci/project/github/poliastro/poliastro/master.svg?style=flat-square&logo=circleci
   :target: https://circleci.com/gh/poliastro/poliastro

.. |travisci| image:: https://img.shields.io/travis/poliastro/poliastro/master.svg?style=flat-square&logo=travis
   :target: https://travis-ci.org/poliastro/poliastro

.. |appveyor| image:: https://img.shields.io/appveyor/ci/Juanlu001/poliastro/master.svg?style=flat-square&logo=appveyor
   :target: https://ci.appveyor.com/project/Juanlu001/poliastro/branch/master

.. |codecov| image:: https://img.shields.io/codecov/c/github/poliastro/poliastro.svg?style=flat-square
   :target: https://codecov.io/github/poliastro/poliastro?branch=master

.. |codeclimate| image:: https://api.codeclimate.com/v1/badges/fd2aa5bf8c4b7984d11b/maintainability
   :target: https://codeclimate.com/github/poliastro/poliastro/maintainability

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://docs.poliastro.space/en/latest/?badge=latest

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/poliastro/poliastro/raw/master/COPYING

.. |doi| image:: https://zenodo.org/badge/11178845.svg?style=flat-square
   :target: https://zenodo.org/badge/latestdoi/11178845

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat-square
   :target: http://www.astropy.org/

.. |mailing| image:: https://img.shields.io/badge/mailing%20list-groups.io-8cbcd1.svg?style=flat-square
   :target: https://groups.io/g/poliastro-dev

.. |matrix| image:: https://img.shields.io/matrix/poliastro:matrix.org.svg?style=flat
   :alt: Join the chat at https://chat.openastronomy.org/#/room/#poliastro:matrix.org
   :target: https://chat.openastronomy.org/#/room/#poliastro:matrix.org

|circleci| |travisci| |appveyor| |codecov| |codeclimate|

|docs| |license| |doi| |astropy| |mailing| |matrix|

poliastro is an open source pure Python package dedicated to problems arising in Astrodynamics and
Orbital Mechanics, such as orbit propagation, solution of the Lambert's
problem, conversion between position and velocity vectors and classical
orbital elements and orbit plotting, focusing on interplanetary applications.
It is released under the MIT license.

.. code-block:: python

    from poliastro.examples import molniya
    from poliastro.plotting import plot

    plot(molniya)

.. image:: https://github.com/poliastro/poliastro/raw/master/docs/source/examples/molniya.png
   :align: center

Documentation
=============

|docs|

Complete documentation, including a user guide and an API reference, can be read on
the wonderful `Read the Docs`_.

https://docs.poliastro.space/

.. _`Read the Docs`: https://readthedocs.org/

Examples
========

.. |mybinder| image:: https://img.shields.io/badge/launch-binder-e66581.svg?style=flat-square
   :target: https://beta.mybinder.org/v2/gh/poliastro/poliastro/master?filepath=index.ipynb


|mybinder|

In the examples directory you can find several Jupyter notebooks with specific
applications of poliastro. You can launch a cloud Jupyter server using `binder`_ to edit
the notebooks without installing anything. Try it out!

https://beta.mybinder.org/v2/gh/poliastro/poliastro/master?filepath=index.ipynb

.. _binder: https://beta.mybinder.org/

Requirements
============

poliastro requires the following Python packages:

* NumPy, for basic numerical routines
* Astropy, for physical units and time handling
* numba (optional), for accelerating the code
* jplephem, for the planetary ephemerides using SPICE kernels
* matplotlib, for orbit plotting
* plotly, for 2D and 3D interactive orbit plotting
* SciPy, for root finding and numerical propagation

poliastro is usually tested on Linux, Windows and OS X on Python
3.5, 3.6 and 3.7 against latest NumPy.

==============  ============  ===================
Platform        Site          Status
==============  ============  ===================
Linux           CircleCI      |circleci|
OS X            Travis CI     |travisci|
Windows x64     Appveyor      |appveyor|
==============  ============  ===================

Installation
============

The easiest and fastest way to get the package up and running is to
install poliastro using `conda <http://conda.io>`_::

  $ conda install poliastro --channel conda-forge

Please check out the `documentation for alternative installation methods`_.

.. _`documentation for alternative installation methods`: https://docs.poliastro.space/en/latest/getting_started.html#alternative-installation-methods

Testing
=======

|codecov|

If installed correctly, the tests can be run using pytest::

  $ python -c "import poliastro.testing; poliastro.testing.test()"
  ===================================== test session starts =====================================
  platform linux -- Python 3.7.1, pytest-4.2.0, py-1.7.0, pluggy-0.8.1
  rootdir: /home/juanlu/.miniconda36/envs/_test37/lib/python3.7/site-packages/poliastro, inifile:
  collected 747 items
  [...]
  ========= 738 passed, 3 skipped, 5 xfailed, 1 xpassed, 13 warnings in 392.12 seconds ==========
  $

Problems
========

If the installation fails or you find something that doesn't work as expected,
please open an issue in the `issue tracker`_.

.. _`issue tracker`: https://github.com/poliastro/poliastro/issues

Contributing
============

.. image:: https://img.shields.io/waffle/label/poliastro/poliastro/1%20-%20Ready.svg?style=flat-square
   :target: https://waffle.io/poliastro/poliastro
   :alt: 'Stories in Ready'

poliastro is a community project, hence all contributions are more than
welcome! For more information, head to `CONTRIBUTING.rst`_.

.. _`CONTRIBUTING.rst`: https://github.com/poliastro/poliastro/blob/master/CONTRIBUTING.rst

Support
=======

|mailing|

Release announcements and general discussion take place on our `mailing list`_.
Feel free to join!

.. _`mailing list`: https://groups.io/g/poliastro-dev

https://groups.io/g/poliastro-dev

Citing
======

If you use poliastro on your project, please
`drop me a line <mailto:juanlu001@gmail.com>`_.

You can also use the DOI to cite it in your publications. This is the latest
one:

|doi|

And this is an example citation format::

 Juan Luis Cano Rodríguez et al.. (2015). poliastro: poliastro 0.4.0. Zenodo. 10.5281/zenodo.17462

License
=======

|license|

poliastro is released under the MIT license, hence allowing commercial
use of the library. Please refer to the COPYING file.

FAQ
===

What's up with the name?
------------------------

poliastro comes from Polimi, which is the shortened name of the Politecnico di
Milano, the Italian university where I was studying while writing this
software. It's my tiny tribute to a place I came to love. *Grazie mille!*

Can I do <insert awesome thing> with poliastro?
-----------------------------------------------

poliastro is focused on interplanetary applications. This has two consequences:

* It tries to be more general than other Flight Dynamics core libraries more
  focused on Earth satellites (see `Related software`_ for a brief list),
  allowing the algorithms to work also for orbits around non-Earth bodies.
* It leaves out certain features that would be too Earth-specific, such as
  TLE reading, SGP4 propagation, groundtrack plotting and others.

.. _`Related software`: https://docs.poliastro.space/en/latest/about.html#related-software

What's the future of the project?
---------------------------------

poliastro is actively maintained and receiving an influx of new contributors
thanks to the generous sponsorship of Google and the European Space Agency.
The best way to get an idea of the roadmap is to see the `Milestones`_ of
the project.

.. _`Milestones`: https://github.com/poliastro/poliastro/milestones
