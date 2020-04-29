Perturbations module
====================

.. graphviz::

   digraph {
      "poliastro.core.perturbations" -> "J2_perturbation", "J3_perturbation", "atmospheric_drag_exponential", "shadow_function", "third_body", "radiation_pressure";
   }

\
\
The :py:mod:`poliastro.core.perturbations` enables the user of :py:mod:`poliastro` to reproduce some
perturbations, since certain conditions can not always be assumed ideal. This module enables to recreate those
non-ideal conditions.

.. automodule:: poliastro.core.perturbations
    :members:
