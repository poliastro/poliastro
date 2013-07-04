=========
poliastro
=========

Overview
--------

These are some Python wrappers to Fortran and MATLAB subroutines useful in
Orbital Mechanics, most of them coming from the `companion software of
Vallado`__.

.. __: http://celestrak.com/software/vallado-sw.asp

Some of them were slightly modified due to errors in the build process,
the handling of relative errors in comparisons and to make them more
suitable to use with different gravitational parameters.

Requirements
------------

The software required by poliastro is:

* Python 3
* Octave
* A Fortran compiler, e.g. gfortran

Also, this Python packages must be present in the system:

* NumPy
* SciPy
* oct2py

poliastro has been tested under

* Linux
* Python 3.3

but there is no reason it shouldn't work under Windows or Mac OS X with
properly configured tools (not tested).

Python 2 compatibility might be accomplished with little syntax changes using
`3to2`_.

.. _3to2: https://pypi.python.org/pypi/3to2

Installation
------------

To install poliastro, just type::

  $ python setup.py install

This might require superuser privileges. To install in a local directory::

  $ python setup.py install --user

Remember to put this directory in the :code:`$PYTHONPATH` if it is not already.

It is recommended that you **never ever use sudo** with distutils, pip,
setuptools and friends in Linux because you might seriously break your
system [1_][2_][3_][4_]. Apart from `per user directories`_, other options
are using `virtualenv`_  or `local installations`_.

.. _1: http://wiki.python.org/moin/CheeseShopTutorial#Distutils_Installation
.. _2: http://stackoverflow.com/questions/4314376/how-can-i-install-a-python-egg-file/4314446#comment4690673_4314446
.. _3: http://workaround.org/easy-install-debian
.. _4: http://matplotlib.1069221.n5.nabble.com/Why-is-pip-not-mentioned-in-the-Installation-Documentation-tp39779p39812.html

.. _`per user directories`: http://stackoverflow.com/a/7143496/554319
.. _`virtualenv`: http://pypi.python.org/pypi/virtualenv
.. _`local installations`: http://stackoverflow.com/a/4325047/554319

Testing
-------

If installed correctly, this should work::

  $ python
  >>> import poliastro
  >>> poliastro.test()
  ...
  OK
  ...
  >>> 

.. TODO: NOT MULTIPROCESSING SAFE, due to oct2py
.. TODO: Ask about libraries and extensions
