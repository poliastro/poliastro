Orbit 
=====

This is probably one of the main important modules of :py:mod:`poliastro` since it
enables the user to create the :py:mod:`poliastro.twobody.orbit.Orbit` objects.

User can create or initialize an orbit in a very different ways but the most common ones are the
following:

.. graphviz::

   digraph {
      "poliastro.twobody.orbit.Orbit" -> "from_vectors", "from_classical", "from_equinoctial", "from_body_ephem";
   }

\
\
Howeverm it is also possible to create 'special' orbits such us :py:mod:`poliastro.twobody.orbit.Orbit.circular`
or even :py:mod:`poliastro.twobody.orbit.Orbit.parabolic`.

Furthermore, orbits can be sampled by :py:mod:`poliastro.twobody.orbit.Orbit.sample` or propagated
in time by :py:mod:`poliastro.twobody.orbit.Orbit.propagate`. But even it is possible to apply a maneuver
to an orbit with the :py:mod:`poliastro.twobody.orbit.Orbit.apply_maneuver`.


.. automodule:: poliastro.twobody.orbit
    :members: