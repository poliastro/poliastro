What's new
==========

New in poliastro 0.4.0
----------------------

.. warning:: This version is not released yet.

This is a new major release, focused on improving stability and code quality.
New angle conversion and modified equinoctial elements functions were added
and an important backwards incompatible change was introduced related to
classical orbital elements.

New features:

* **Angle conversion functions**: Finally brought back from poliastro 0.1,
  new functions were added to convert between true \\(\\nu\\), eccentric
  \\(E\\) and mean \\(M\\) anomaly, see `#45`_.
* **Equinoctial elements**: Now it's possible to convert between classical
  and equinoctial elements, as well as from/to position and velocity vectors,
  see `#61`_.
* **Numerical propagation**: A new propagator using
  SciPy Dormand & Prince 8(5,3) integrator was added, see `#64`_.

.. _`#45`: https://github.com/poliastro/poliastro/pull/45
.. _`#61`: https://github.com/poliastro/poliastro/pull/61
.. _`#64`: https://github.com/poliastro/poliastro/pull/64

Other highlights:

* **MIT license**: The project has been relicensed to a more popular license.
  poliastro remains commercial-friendly through a permissive, OSI-approved
  license.
* **Python 3.5 and NumPy 1.10 compatibility**. poliastro retains compatibility
  with legacy Python (Python 2) and NumPy 1.9. *Next version will be Python 3
  only*.

Backward incompatible changes:

* **Creation of orbits from classical elements**: poliastro has
  switched to the *semilatus rectum* \\(p\\) instead of the semimajor
  axis \\(a\\) to define ``State`` objects, and the function has been renamed
  to :py:meth:`~poliastro.twobody.State.from_classical`. Please update your
  programs accordingly.
* Removed specific angular momentum \\(h\\) property to avoid a name clash
  with the fourth modified equinoctial element, use ``norm(ss.h_vec)``
  instead.

New in poliastro 0.3.1
----------------------

This is a new minor release, with some bug fixes backported from the main
development branch.

Bugs fixed:

* Fixed installation problem in Python 2.
* `Issue #49`_: Fix velocity units in ``ephem``.
* `Issue #50`_: Fixed ``ZeroDivisionError`` when propagating with time zero.

.. _`Issue #49`: https://github.com/poliastro/poliastro/issues/49
.. _`Issue #50`: https://github.com/poliastro/poliastro/issues/50

New in poliastro 0.3.0
----------------------

This is a new major release, focused on switching to a pure Python codebase.
Lambert problem solving and ephemerides computation came back, and a couple
of bugs were fixed.

New features:

* **Pure Python codebase**: Forget about Fortran linking problems and
  nightmares on Windows, because now poliastro is a pure Python package.
  A new dependency, numba, was introduced to accelerate the algorithms,
  but poliastro will use it only if it is installed.
* **Lambert problem solving**: New module :py:mod:`~poliastro.iod` to
  determine an orbit given two position vectors and the time of flight.
* `PR #42`_: **Planetary ephemerides computation**: New module
  :py:mod:`~poliastro.ephem` with functions to deal with SPK files and
  compute position and velocity vectors of the planets.
* `PR #38`_: New method :py:meth:`~poliastro.twobody.State.parabolic` to create parabolic orbits.
* New conda package: visit `poliastro binstar channel`_!
* New organization and logo.

.. _`PR #42`: https://github.com/poliastro/poliastro/pull/42
.. _`PR #38`: https://github.com/poliastro/poliastro/pull/38

.. _`poliastro binstar channel`: https://binstar.org/poliastro

Bugs fixed:

* `Issue #19`_: Fixed plotting region for parabolic orbits.
* `Issue #37`_: Fixed creation of parabolic orbits.

.. _`Issue #19`: https://github.com/poliastro/poliastro/issues/19
.. _`Issue #37`: https://github.com/poliastro/poliastro/issues/37

New in poliastro 0.2.1
----------------------

This is a bugfix release, no new features were introduced since 0.2.0.

* Fixed `#35`_ (failing tests with recent astropy versions), thanks to
  Sam Dupree for the bug report.
* Updated for recent Sphinx versions.

.. _`#35`: https://github.com/poliastro/poliastro/issues/35

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

.. _`PR #14`: https://github.com/poliastro/poliastro/pull/14
.. _wiki: https://github.com/poliastro/poliastro/wiki
.. _Plyades: https://github.com/helgee/Plyades

These features were removed temporarily not to block the release and will
see the light again in poliastro 0.3:

* Conversion between anomalies.
* Ephemerides calculations, will look into Skyfield and the JPL ephemerides
  prepared by Brandon Rhodes (see `issue #4`_).
* Lambert problem solving.
* Perturbation analysis.

.. _`issue #4`: https://github.com/poliastro/poliastro/issues/4
