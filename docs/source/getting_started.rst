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
* pytest, for running the tests from the package

poliastro is usually tested on Linux, Windows and OS X on Python
3.5 and 3.6 against latest NumPy.

Installation
------------

The easiest and fastest way to get the package up and running is to
install poliastro using `conda <https://conda.io/docs/>`_::

  $ conda install poliastro --channel conda-forge

.. note::

    We encourage users to use conda and the
    `conda-forge <https://conda-forge.org/>`_ packages for convenience,
    especially when developing on Windows.

If the installation fails for any reason, please open an issue in the
`issue tracker`_.

Alternative installation methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't want to use conda you can `install poliastro from PyPI`_
using pip::

  $ pip install numpy  # Run this one first for pip 9 and older!
  $ pip install poliastro

Finally, you can also install the latest development version of poliastro
`directly from GitHub`_::

  $ pip install https://github.com/poliastro/poliastro/archive/master.zip

This is useful if there is some feature that you want to try, but we did not
release it yet as a stable version. Although you might find some unpolished
details, these development installations should work without problems. If
you find any, please open an issue in the `issue tracker`_.

.. _`install poliastro from PyPI`: https://pypi.python.org/pypi/poliastro/
.. _`directly from GitHub`: http://github.com/poliastro/poliastro

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

If installed correctly, the tests can be run using pytest::

  $ python -c "import poliastro.testing; poliastro.testing.test()"
  Running unit tests for poliastro
  [...]
  OK
  $ 

If for some reason any test fails, please report it in the `issue tracker`_.

.. _`issue tracker`: https://github.com/poliastro/poliastro/issues
