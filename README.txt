=========
poliastro
=========

These are some Python wrappers to Fortran and MATLAB subroutines useful in
Orbital Mechanics, most of them coming from the `companion software of
Vallado`__.

.. __: http://celestrak.com/software/vallado-sw.asp

Some of them were slightly modified due to errors in the build process,
the handling of relative errors in comparisons and to make them more
suitable to use with different gravitational parameters.

The software required by poliastro is:

* Python 3
* Octave
* A Fortran compiler, e.g. gfortran

Also, this Python packages must be present in the system:

* NumPy
* SciPy
* oct2py

To install poliastro, just type::

  $ python setup.py install

This might require superuser privileges. To install in a local directory::

  $ python setup.py install --prefix=~/.local

Remember to put this directory in the :code:`$PYTHONPATH` if it is not already.

If installed correctly, this should work::

  $ python
  >>> import poliastro
  >>> dir(poliastro)
  ['M2nu', '__builtins__', '__cached__', ...
  'rv2coe', 'target', 'twobody', 'util']

poliastro has been tested under

* Linux
* Python 3.3

but there is no reason for it not to work under Windows or Mac OS X.
Python 2 support might be accomplished with little syntax changes using
`3to2`_.

.. _3to2: https://pypi.python.org/pypi/3to2

.. TODO: Change argument order
.. TODO: Test in different systems
.. TODO: Add license
.. TODO: Release

.. TODO: Ask about libraries and extensions
