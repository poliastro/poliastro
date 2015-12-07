Getting started
===============

Requirements
------------

poliastro requires the following Python packages:

* NumPy, for basic numerical routines
* Astropy, for physical units and time handling
* numba (optional), for accelerating the code
* jplephem, for the planetary ephemerides using SPICE kernels
* matplotlib, for orbit plotting
* scipy, for root finding and numerical propagation

poliastro is usually tested on Linux, Windows and OS X on Python 2.7, 3.3,
3.4 and 3.5, using NumPy 1.9 and 1.10 (single codebase).

Installation
------------

The easiest and fastest way to get the package up and running is to
install poliastro using `conda <http://conda.io>`_::

  $ conda install poliastro --channel poliastro

You can also `install poliastro from PyPI`_ using pip, given that you already
have all the requirements::

  $ pip install poliastro

You can also `download poliastro source from GitHub`_ and type::

  $ python setup.py install

Development installations are also supported thanks to setuptools::

  $ python setup.py develop

.. _`install poliastro from PyPI`: https://pypi.python.org/pypi/poliastro/
.. _`download poliastro source from GitHub`: http://github.com/poliastro/poliastro

.. warning::

    It is recommended that you **never ever use sudo** with distutils, pip,
    setuptools and friends in Linux because you might seriously break your
    system [1_][2_][3_][4_]. Options are `per user directories`_, `virtualenv`_
    or `local installations`_.

.. _1: http://wiki.python.org/moin/CheeseShopTutorial#Distutils_Installation
.. _2: http://stackoverflow.com/questions/4314376/how-can-i-install-a-python-egg-file/4314446#comment4690673_4314446
.. _3: http://workaround.org/easy-install-debian
.. _4: http://matplotlib.1069221.n5.nabble.com/Why-is-pip-not-mentioned-in-the-Installation-Documentation-tp39779p39812.html

.. _`per user directories`: http://stackoverflow.com/a/7143496/554319
.. _`virtualenv`: http://pypi.python.org/pypi/virtualenv
.. _`local installations`: http://stackoverflow.com/a/4325047/554319

Testing
-------

If installed correctly, the tests can be run using py.test::

  $ python -c "import poliastro; poliastro.test()"
  Running unit tests for poliastro
  [...]
  OK
  $ 

If for some reason any test fails, please report it in the `issue tracker`_.

.. _`issue tracker`: https://github.com/poliastro/poliastro/issues
