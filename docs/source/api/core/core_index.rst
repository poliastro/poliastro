Core level API
==============

The :py:mod:`poliastro.core` includes a set of modules that form the kernel of :py:mod:`poliastro`. All defined functions
in these modules work with `Numba’s JIT compiler` (`Numba web site <http://numba.pydata.org/>`__), since many of the computations
behind poliastro involve hard numerical methods. This is the Low-Level API of :py:mod:`poliastro`, being the basis of the High-Level API,
the one expected to be used.

.. graphviz::

   digraph {
      "poliastro.core" -> "angles", "elements", "hyper", "iod", "perturbations", "propagation", "stumpff", "util";
   }

.. toctree::
    :hidden:
    :maxdepth: 1

    angles
    elements
    hyper
    iod
    perturbations
    propagation
    stumpff
    util

