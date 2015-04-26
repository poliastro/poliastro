What's new
==========

New in poliastro 0.2.1
----------------------

This is a bugfix release, no new features were introduced since 0.2.0.

* Fixed `#35`_ (failing tests with recent astropy versions), thanks to
  Sam Dupree for the bug report.
* Updated for recent Sphinx versions.

.. _`#35`: https://github.com/Pybonacci/poliastro/issues/35

New in poliastro 0.2
--------------------

* **Totally refactored code** to provide a more pythonic API (see `PR #14`_
  and `wiki`_ for further information) heavily inspired by `Plyades`_ by
  Helge Eichhorn.

  * Mandatory use of **physical units** through :code:`astropy.units`.
  * Object-oriented approach: :py:class:`~poliastro.twobody.State` and
    :py:class:`~poliastro.maneuver.Maneuver` classes.
  * Vector quantities: results not only have magnitude now, but also direction
    (see for example maneuvers).

* Easy plotting of orbits in two dimensions using matplotlib.
* Module :code:`example` with sample data to start testing the library.

.. _`PR #14`: https://github.com/Pybonacci/poliastro/pull/14
.. _wiki: https://github.com/Pybonacci/poliastro/wiki
.. _Plyades: https://github.com/helgee/Plyades

These features were removed temporarily not to block the release and will
see the light again in poliastro 0.3:

* Conversion between anomalies.
* Ephemerides calculations, will look into Skyfield and the JPL ephemerides
  prepared by Brandon Rhodes (see `issue #4`_).
* Lambert problem solving.
* Perturbation analysis.

.. _`issue #4`: https://github.com/Pybonacci/poliastro/issues/4
