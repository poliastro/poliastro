Initial Orbit Determination (IOD) module
========================================

The initial orbit determination problem consists on finding the orbit
that passes trough two different position vectors during a given finite
period of time.

:py:mod:`poliastro` includes two algorithms under the name 'lambert' for solving this problem, also known as
Lambert's problem. :py:mod:`poliastro.iod` is composed by the following two sub-modules:

.. toctree::
    :maxdepth: 2

    izzo
    vallado