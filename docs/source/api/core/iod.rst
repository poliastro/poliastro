Initial Orbit Determination (IOD) module
========================================

.. graphviz::

   digraph {
      "poliastro.core.iod" -> "vallado", "izzo";
   }

\
\
The :py:mod:`polaistro.core.iod` deals with the problem of determining an orbit
being given two position vectors along it and the that the body takes to travel from
one to another. 

This problem is known as Lambert's problem and many algorithms have developed to solved since
the main difficult of it is focused on numerical methods. This module contains two different
functions that enable to solve the Lambert problem.


.. automodule:: poliastro.core.iod
   :members:
