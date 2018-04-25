What's new
==========

poliastro 0.9.0 - 2018-04-25
----------------------------

This major release received lots of improvements in the 2D plotting code
and propagation functions, introduced the new perturbation framework and
paved the way for the `Python in Astronomy 2018`_ workshop and the
`Google Summer of Code 2018`_ program.

.. _`Python in Astronomy 2018`: http://openastronomy.org/pyastro/2018/
.. _`Google Summer of Code 2018`: https://summerofcode.withgoogle.com/

New features
............

* **New experimental 2D Plotly backend**: A new :py:class:`~poliastro.plotting.OrbitPlotter2D`
  class was introduced that uses Plotly instead of matplotlib for the rendering. There are
  still some issues that should be resolved when we take advantage of the latest Plotly version,
  hence the "experimental" nature.
* **New propagators**: A new Keplerian propagator :py:meth:`~poliastro.twobody.propagation.mean_motion`
  was introduced that has better convergence properties than :py:meth:`~poliastro.twobody.propagation.kepler`,
  so now the user can choose.
* **New perturbation functions**: A new module :py:mod:`poliastro.twobody.perturbations` was introduced
  that contains perturbation accelerations that can be readily used with
  :py:meth:`~poliastro.twobody.propagation.cowell`. So far we implemented J2 and atmospheric drag effects,
  and we will add more during the summer. Check out the User Guide for examples!
* **Support for different propagators in sampling**: With the introduction of new propagators and perturbation
  accelerations, now the user can easily sample over a period of time using any of them. We are eager to see
  what experiments you come up with!
* **Easy plotting of the Solar System**: A new function :py:meth:`~poliastro.plotting.plot_solar_system` was
  added to easily visualize our inner or complete Solar System in 2D plots.

Other highlights
................

* **poliastro participates in Google Summer of Code thanks to OpenAstronomy!**
  More information `in the poliastro blog <http://blog.poliastro.space/2018/02/22/2018-02-22-join-poliastro-google-summer-of-code/>`_.
* **poliastro will be presented at the Python in Astronomy 2018 workshop** to be held
  at Center for Computational Astrophysics at the Flatiron Institute in New York, USA.
  You can read `more details about the event here <http://openastronomy.org/pyastro/2018/>`_.

New contributors
................

This is the complete list of the people that contributed to this release,
with a + sign indicating first contribution.

* Juan Luis Cano
* Pablo Galindo+
* Matt Ettus+
* Shreyas Bapat+
* Ritiek Malhotra+
* Nikita Astrakhantsev+

Bugs fixed:
...........

* `Issue #294`_: Default steps 2D plots were too visible

.. _`Issue #294`: https://github.com/poliastro/poliastro/issues/294

Backward incompatible changes
.............................

* Now the :py:meth:`poliastro.twobody.Orbit.sample` method returns a tuple of (times, positions).
* All the propagator methods changed their signature
  and now accept :py:class:`~poliastro.twobody.Orbit` objects.

poliastro 0.8.0 - 2017-11-18
----------------------------

This is a new major release, focused on bringing 3D plotting functions and
preparing the material for the Open Source Cubesat Workshop.

New features
............

* **Sampling method** for :py:class:`~poliastro.twobody.Orbit` objects that returns
  an array of positions. This was already done in the plotting functions and will
  help providing other applications, such as exporting an Orbit to other formats.
* **3D plotting functions**: finally poliastro features a new high level object,
  :py:class:`poliastro.plotting.OrbitPlotter3D`, that uses Plotly to represent
  orbit and trajectories in 3D. The venerable notebook about the trajectory of
  rover Curiosity has been updated accordingly.
* **Propagation to a certain date**: now apart from specifying the total elapsed
  time for propagation or time of flight, we can directly specify a target date
  in :py:meth:`poliastro.twobody.orbit.Orbit.propagate`.
* **Hyperbolic anomaly conversion**: we implemented the conversion of hyperbolic
  to mean and true anomaly to complement the existing eccentric anomaly functions
  and improve the handling of hyperbolic orbits in :py:mod:`poliastro.twobody.angles`.

Other highlights
................

* **poliastro is now an Astropy affiliated package**, which gives the project a privileged
  position in the Python ecosystem. Thank you, Astropy core developers! You can read
  `the evaluation here <https://github.com/poliastro/poliastro/issues/279>`_.
* **poliastro will be presented at the first Open Source Cubesat Workshop** to be held
  at the European Space Operations Centre in Darmstadt, Germany. You can read
  `the full program of the event here <http://oscw.space/>`_.

New contributors
................

This is the complete list of the people that contributed to this release,
with a + sign indicating first contribution.

* Juan Luis Cano
* Antonio Hidalgo
* mattrossman+
* Roshan Jossey+

Bugs fixed:
...........

* `Issue #275`_: Converting from true to mean anomaly fails for hyperbolic orbits

.. _`Issue #275`: https://github.com/poliastro/poliastro/issues/275

Backward incompatible changes
.............................

* The :code:`ephem` module has been removed in favor of the
  :code:`astropy.coordinates.get_body_barycentric_posvel` function.

poliastro 0.7.0 - 2017-09-15
----------------------------

This is a new major release, which adds new packages and modules,
besides fixing several issues.

New features:
.............

* **NEOS package**: a new package has been added to poliastro, :py:mod:`~poliastro.neos`
  package. It provides several ways of getting NEOs (Near Earth Objects) data from NASA
  databases, online and offline.
* **New patched conics module**. New module containing a function to compute
  the radius of the Sphere of Influence (SOI).
* **Use Astropy for body ephemerides**. Instead of downloading the SPK
  files ourselves, now we use Astropy builtin capabilities. This also
  allows the user to select a builtin ephemerides that does not require
  external downloads. See `#131`_ for details.
* **Coordinates and frames modules**: new modules containing transformations between ICRS
  and body-centered frame, and perifocal to body_centered, :py:mod:`~poliastro.coordinates`
  as well as Heliocentric coordinate frame in :py:mod:`~poliastro.frames` based on Astropy
  for NEOs.
* **Pip packaging**: troublesome dependencies have been released in wheel format,
  so poliastro can now be installed using pip from all platforms.
* **Legend plotting**: now label and epoch are in a figure legend, which ends with
  the ambiguity of the epochs when having several plots in the same figure.

.. _`#131`: https://github.com/poliastro/poliastro/issues/131


Other highlights:
.................

* **Joined Open Astronomy**: we are now part of `Open Astronomy`_, a
  collaboration between open source astronomy and astrophysics projects
  to share resources, ideas, and to improve code.
* **New constants module**: poliastro has now a :py:mod:`~poliastro.constants` module,
  with GMs and radii of solar system bodies.
* **Added Jupyter examples**: poliastro examples are now available in the
  documentation as Jupyter notebooks, thanks to `nbsphinx`_.
* **New Code of Conduct**: poliastro community now has a Code of conduct.
* **Documentation update**: documentation has been updated with new installation
  ways, propagation and NEOs examples, "refactored" code and images, improved contribution
  guidelines and intersphinx extension.
* **New success stories**: two new success stories have been added to documentation.
* **Bodies now have a parent**. It is now possible to specify the attractor
  of a body.
* **Relative definition of Bodies**. Now it is possible to define Body parameters
  with respect to another body, and also add any number of properties in a simple
  way.

.. _`nbsphinx`: http://nbsphinx.readthedocs.io/en/latest/
.. _`Open Astronomy`: http://openastronomy.org/members/

New contributors
................

Thanks to the generous SOCIS grant from the European Space Agency,
Antonio Hidalgo has devoted three months developing poliastro full time
and gained write acces to the repository.

This is the complete list of the people that contributed to this release,
with a + sign indicating first contribution.

* Juan Luis Cano
* MiguelHB+
* Antonio Hidalgo+
* Zac Miller+
* Fran Navarro+
* Pablo Rodríguez Robles+

Bugs fixed:
...........

* `Issue #205`_: Bug when plotting orbits with different epochs.
* `Issue #128`_: Missing ephemerides if no files on import time.
* `Issue #131`_: Slightly incorrect ephemerides results due to improper time scale.
* `Issue #130`_: Wrong attractor size when plotting different orbits.

.. _`Issue #205`: https://github.com/poliastro/poliastro/issues/205
.. _`Issue #128`: https://github.com/poliastro/poliastro/issues/128
.. _`Issue #131`: https://github.com/poliastro/poliastro/issues/131
.. _`Issue #130`: https://github.com/poliastro/poliastro/issues/130

Backward incompatible changes:
..............................

* **Non-osculating orbits**: removed support for non-osculating orbits.
  :code:`plotting.plot()` calls containing :code:`osculating` parameter should be
  replaced.

poliastro 0.6.0 - 2017-02-12
----------------------------

This major release was focused on refactoring some internal core
parts and improving the propagation functionality.

Highlights:
...........

* **Support Python 3.6**. See `#144`_.
* **Introduced ``Orbit`` objects** to replace ``State`` ones. The latter
  has been simplified, reducing some functionality, now their API
  has been moved to the former. See the User Guide and the examples for
  updated explanations. See `#135`_.
* **Allow propagation functions to receive a callback**. This paves the
  way for better plotting and storage of results. See `#140`_.

.. _`#135`: https://github.com/poliastro/poliastro/pull/135
.. _`#140`: https://github.com/poliastro/poliastro/pull/140
.. _`#144`: https://github.com/poliastro/poliastro/pull/144

poliastro 0.5.0 - 2016-03-06
----------------------------

This is a new major release, focused on expanding the initial orbit
determination capabilities and solving some infrastructure challenges.

New features:
.............

* **Izzo's algorithm for the Lambert problem**: Thanks to this algorithm
  multirevolution solutions are also returned. The old algorithm is kept
  on a separate module.

Other highlights:
.................

* **Documentation on Read the Docs**: You can now browse previous releases
  of the package and easily switch between released and development versions.
* **Mailing list**: poliastro now has a mailing list hosted on groups.io.
  Come and join!
* **Clarified scope**: poliastro will now be focused on interplanetary
  applications, leaving other features to the new `python-astrodynamics`_
  project.

.. _`python-astrodynamics`: http://python-astrodynamics.org/

Bugs fixed:
...........

* `Issue #110`_: Bug when plotting State with non canonical units

.. _`Issue #110`: https://github.com/poliastro/poliastro/issues/110

Backward incompatible changes:
..............................

* **Drop Legacy Python**: poliastro 0.5.x and later will support only
  Python 3.x. We recommend our potential users to create dedicated virtual
  environments using conda or virtualenv or to contact the developers to fund
  Python 2 support.
* **Change ``lambert`` function API**: The functions for solving Lambert's
  problem are now _generators_, even in the single revolution case.
  Check out the User Guide for specific examples.
* **Creation of orbits from classical elements**: poliastro has
  reverted the switch to the *semilatus rectum* \\(p\\) instead of the semimajor
  axis \\(a\\) made in 0.4.0, so \\(a\\) must be used again. This change is
  definitive.

poliastro 0.4.2 - 2015-12-24
----------------------------

Fixed packaging problems.

poliastro 0.4.0 - 2015-12-13
----------------------------

This is a new major release, focused on improving stability and code quality.
New angle conversion and modified equinoctial elements functions were added
and an important backwards incompatible change was introduced related to
classical orbital elements.

New features:
.............

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
.................

* **MIT license**: The project has been relicensed to a more popular license.
  poliastro remains commercial-friendly through a permissive, OSI-approved
  license.
* **Python 3.5 and NumPy 1.10 compatibility**. poliastro retains compatibility
  with legacy Python (Python 2) and NumPy 1.9. *Next version will be Python 3
  only*.

Bugs fixed:
...........

* `Issue #62`_: Conversion between coe and rv is not transitive
* `Issue #69`_: Incorrect plotting of certain closed orbits

.. _`Issue #62`: https://github.com/poliastro/poliastro/issues/62
.. _`Issue #69`: https://github.com/poliastro/poliastro/issues/69

Backward incompatible changes:
..............................

* **Creation of orbits from classical elements**: poliastro has
  switched to the *semilatus rectum* \\(p\\) instead of the semimajor
  axis \\(a\\) to define ``State`` objects, and the function has been renamed
  to :py:meth:`~poliastro.twobody.State.from_classical`. Please update your
  programs accordingly.
* Removed specific angular momentum \\(h\\) property to avoid a name clash
  with the fourth modified equinoctial element, use ``norm(ss.h_vec)``
  instead.

poliastro 0.3.1 - 2015-06-30
----------------------------

This is a new minor release, with some bug fixes backported from the main
development branch.

Bugs fixed:
...........

* Fixed installation problem in Python 2.
* `Issue #49`_: Fix velocity units in ``ephem``.
* `Issue #50`_: Fixed ``ZeroDivisionError`` when propagating with time zero.

.. _`Issue #49`: https://github.com/poliastro/poliastro/issues/49
.. _`Issue #50`: https://github.com/poliastro/poliastro/issues/50

poliastro 0.3.0 - 2015-05-09
----------------------------

This is a new major release, focused on switching to a pure Python codebase.
Lambert problem solving and ephemerides computation came back, and a couple
of bugs were fixed.

New features:
.............

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
...........

* `Issue #19`_: Fixed plotting region for parabolic orbits.
* `Issue #37`_: Fixed creation of parabolic orbits.

.. _`Issue #19`: https://github.com/poliastro/poliastro/issues/19
.. _`Issue #37`: https://github.com/poliastro/poliastro/issues/37

poliastro 0.2.1 - 2015-04-26
----------------------------

This is a bugfix release, no new features were introduced since 0.2.0.

* Fixed `#35`_ (failing tests with recent astropy versions), thanks to
  Sam Dupree for the bug report.
* Updated for recent Sphinx versions.

.. _`#35`: https://github.com/poliastro/poliastro/issues/35

poliastro 0.2 - 2014-08-16
--------------------------

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
