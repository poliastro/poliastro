Twobody module
==============

The :py:mod:`poliastro.twobody` is one of the most important modules since it includes the
:py:mod:`polastro.twobody`

Submodules
----------

.. toctree::
    :maxdepth: 1

    angles
    decorators
    elements
    events
    mean_elements
    propagation
    sampling

Orbit
-----

Users can create or initialize an orbit in a very different ways but the most common ones are the
following:

.. graphviz::

   digraph {
      "poliastro.twobody.orbit.Orbit" -> "from_vectors", "from_classical", "from_equinoctial", "from_ephem";
   }

\
\
It is also possible to create 'special' orbits such us :py:mod:`poliastro.twobody.orbit.Orbit.circular`
or even :py:mod:`poliastro.twobody.orbit.Orbit.parabolic`.

Furthermore, orbits can be sampled by :py:mod:`poliastro.twobody.orbit.Orbit.sample` or propagated
in time by :py:mod:`poliastro.twobody.orbit.Orbit.propagate`. It is also possible to apply a maneuver
to an orbit with the :py:mod:`poliastro.twobody.orbit.Orbit.apply_maneuver`.

.. automodule:: poliastro.twobody.orbit
    :members:
