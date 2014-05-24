Getting started
===============

Requirements
------------

poliastro requires NumPy, Astropy and matplotlib and is usually tested on Linux
on Python 2.7 and Python 3.3 (single codebase).

Besides, a Fortran compiler is needed to build the extensions. poliastro
is usually built with gfortran.

There is no reason it shouldn't work under Windows or Mac OS X with
properly configured tools (not tested).

Installation
------------

To install poliastro just use pip::

  $ pip install poliastro

To install poliastro from source, just type::

  $ python setup.py install

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

