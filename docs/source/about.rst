About poliastro
===============

Overview
--------

**poliastro** is an open source collection of Python subroutines for solving
problems in Astrodynamics and Orbital Mechanics.

poliastro combines cutting edge technologies like Python JIT compiling
(using numba) with young, well developed astronomy packages (like astropy and
jplephem) to provide a user friendly API for solving Astrodynamics problems.
It is therefore a experiment to mix the best Python open source practices
with my love for Orbital Mechanics.

Since I have only solved easy academic problems I cannot assess the
suitability of the library for professional environments, though I am aware
that at least a company that uses it.

History
-------

I started poliastro as a wrapper of some MATLAB and Fortran algorithms that I
needed for a University project: having good performance was a must, so pure
Python was not an option. As a three language project, it was only known to
work in my computer, and I had to fight against oct2py and f2py for long
hours.

Later on, I enhanced poliastro plotting capabilities to serve me in further
University tasks. I removed the MATLAB (Octave) code and kept only the
Fortran algorithms. Finally, when numba was mature enough, I implemented
everything in pure Python and poliastro 0.3 was born.

Related software
----------------

These are some projects which share similarities with poliastro or which
served as inspiration:

* `astropy`_: According to its website, "The Astropy Project is a community
  effort to develop a single core package for Astronomy in Python and foster
  interoperability between Python astronomy packages". Not only it provides
  important core features for poliastro like time and physical units handling,
  but also sets a high bar for code quality and documentation standards. A
  truly inspiring project.
* `Skyfield`_: Another Astronomy Python package focused on computing
  observations of planetary bodies and Earth satellites written by Brandon
  Rhodes. It is the successor of pyephem, also written by him, but skyfield
  is a pure Python package and provides a much cleaner API.
* `Plyades`_: A pioneering astrodynamics library written in Python by Helgee
  Eichhorn. Its clean and user friendly API inspired me to completely refactor
  poliastro 0.2 so it could be much easier to use. It has been stalled for
  a while, but at the moment of writing these lines its author is pushing new
  commits.
* `orbital`_: Yet another orbital mechanics Python library written by Frazer
  McLean. It is very similar to poliastro (orbital plotting module was
  inspired in mine) but its internal structure is way smarter. It is more
  focused in plotting and it even provides 3D plots and animations.
* `orekit-python-wrapper`_: According to its website, "The Orekit python
  wrapper enables to use Orekit within a normal python environment", using
  JCC. Orekit is a well-stablished, mature open source library for
  Astrodynamics written in Java strongly supported by several space agencies.
  The Python wrapper is developed by the Swedish Space Corporation.

.. _astropy: http://www.astropy.org/
.. _Skyfield: http://rhodesmill.org/skyfield/
.. _Plyades: http://plyades.readthedocs.org/en/latest/
.. _orbital: http://pythonhosted.org/OrbitalPy/
.. _orekit-python-wrapper: https://www.orekit.org/forge/projects/orekit-python-wrapper/wiki

Future ideas
------------

These are some things that I would love to implement in poliastro to expand
its capabilities:

* 3D plotting of orbits
* Continuous thrust maneuvers
* Tisserand graphs
* Porkchop plots

Note of the original author
---------------------------

I am Juan Luis Cano Rodríguez (two names and two surnames, it's the Spanish
way!), an Aerospace Engineer with a passion for Astrodynamics
and the Open Source world. Before poliastro started to be a truly community
project, I started it when I was an Erasmus student
at Politecnico di Milano, an important technical university in Italy which
deeply influenced my life and ambitions and gave name to the library itself.
It is and always will be my tiny tribute to a country that will always be in
my heart and to people that never ceased to inspire me. *Grazie mille!*
