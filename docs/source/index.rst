poliastro - Astrodynamics in Python
===================================

.. figure:: _static/molniya.png
   :align: right
   :figwidth: 300
   :alt: Molniya orbit

   Plot of a Molniya orbit around the Earth
   (\\(a = 26600\\,\\mathrm{km}, e = 0.75,
   i = 63.4 \\mathrm{{}^{\\circ}} \\)).

poliastro is a collection of Python wrappers to Fortran subroutines useful in
Orbital Mechanics, such as:

* Orbit propagation
* Conversion between position and velocity vectors and classical orbital
  elements
* Hohmann and bielliptic maneuvers computation
* Orbit plotting

And more to come!

.. code-block:: python

    from poliastro.examples import molniya
    
    molniya.plot()

Contents
--------

.. toctree::
   :maxdepth: 2

   getting_started
   api

