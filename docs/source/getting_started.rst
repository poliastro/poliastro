Getting started
===============

Requirements
------------

poliastro requires a number of Python packages, notably:

* Astropy, for physical units and time handling
* NumPy, for basic numerical routines
* jplephem, for the planetary ephemerides using SPICE kernels
* matplotlib, for static orbit plotting
* numba (when using CPython), for accelerating the code
* Plotly, for interactive orbit plotting
* SciPy, for root finding and numerical propagation

poliastro is usually tested on Linux and Windows on Python
3.7 and 3.8 against latest NumPy.
It should work on OS X without problems.

Installation
------------

The easiest and fastest way to get the package up and running is to
install poliastro using `conda <https://conda.io/docs/>`_::

  $ conda install -c conda-forge poliastro=0.14

.. note::

    We encourage users to use conda and the
    `conda-forge <https://conda-forge.org/>`_ packages for convenience,
    especially when developing on Windows.
    It is recommended to create a new environment.

If the installation fails for any reason, please open an issue in the
`issue tracker`_.

Alternative installation methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you don't want to use conda you can `install poliastro from PyPI`_
using pip::

  $ pip install numpy  # Run this one first for pip 9 and older!
  $ pip install poliastro[jupyter] pytest

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

Using poliastro on JupyterLab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the release of Plotly 3.0, plotting orbits using poliastro is easier than ever.

You have to install three extensions of JupyterLab to make your experience smooth::

  $ jupyter labextension install @jupyter-widgets/jupyterlab-manager@0.38
  $ jupyter labextension install plotlywidget@0.7.0
  $ jupyter labextension install @jupyterlab/plotly-extension@0.18.1

And as the documentation of JupyterLab Extensions states:

  "In order to install JupyterLab extensions, you need to have Node.js version 4 or later installed."

If you face any further issues, you can refer to the `installation guide by Plotly`_.

.. _`installation guide by Plotly`: https://github.com/plotly/plotly.py/blob/master/README.md#jupyterlab-support-python-35

Problems and suggestions
------------------------

If for any reason you get an unexpected error message or an incorrect result,
or you want to let the developers know about your use case,
please open a new issue in the `issue tracker`_ and we will try to answer promptly.

.. _`issue tracker`: https://github.com/poliastro/poliastro/issues
