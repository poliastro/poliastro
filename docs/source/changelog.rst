What's new
==========

New in poliastro 0.2
--------------------

* **Totally refactored code** to provide a more pythonic API (see `PR #14`_
  and `wiki`_ for further information).
* Mandatory use of **physical units** through :code:`astropy.units`.
* Easy plotting of orbits in two dimensions using matplotlib.

.. _`PR #14`: https://github.com/Pybonacci/poliastro/pull/14
.. _wiki: https://github.com/Pybonacci/poliastro/wiki

These features were removed temporarily not to block the release and will
see the light again in poliastro 0.3:

* Conversion between anomalies.
* Ephemerides calculations, will look into Skyfield and the JPL ephemerides
  prepared by Brandon Rhodes (see `issue #4`_).
* Lambert problem solving.
* Perturbation analysis.

.. _`issue #4`: https://github.com/Pybonacci/poliastro/issues/4
