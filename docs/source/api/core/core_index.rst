Core level API
==============

The :guilabel:`poliastro.core` includes a set of modules that form the kernel of :guilabel:`poliastro`. All defined functions
in these modules work with :guilabel:`Numba’s JIT compiler` (`Numba web site <http://numba.pydata.org/>`__), since many of the computations
behind poliastro involve hard numerical methods.

This is the  :guilabel:`Low-Level API` of :guilabel:`poliastro`, being the basis of the :guilabel:`High-Level API`, the one expected to be
used. The modules included are listed below:

.. toctree::
    :maxdepth: 1

    angles
    elements
    hyper
    iod
    perturbations
    propagators
    stumpff
    util

