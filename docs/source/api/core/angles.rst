Angles module
=============

.. graphviz::

   digraph {
      "poliastro.core.angles" -> "D_to_nu", "nu_to_D", "nu_to_E", "nu_to_F", "E_to_nu", "F_to_nu", 
                                 "M_to_E", "M_to_F", "M_to_D", "E_to_M", "F_to_M", "D_to_M", "M_to_nu",
                                 "nu_to_M", "fp_angle";
   }

\
\
The :py:mod:`poliastro.core.angles` module contains functions related to conversion between different angles used
to define different orbital elements.

.. automodule:: poliastro.core.angles
    :members:
