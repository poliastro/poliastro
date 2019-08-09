High level API
==============

The :py:mod:`poliastro.twobody` is on of the main modules of :py:mod:`poliastro` since enables
the user in a very intuitive way to create, propagate and apply maneuvers to an orbit due to its
OOP nature.


.. graphviz:: 
    
    digraph {
        "poliastro" -> "twobody", "threebody", "bodies", "neos", "plotting", "iod", "constants", "coordinates", "examples", "frames",
                       "maneuver"

    }

.. toctree::
    :hidden:
    :maxdepth: 2

    atmosphere/atmosphere_index
    twobody/twobody_index
    threebody/threebody_index
    bodies
    neos/neos_index
    plotting/plotting_index
    iod/iod_index
    constants
    coordinates
    examples
    frames
    maneuver
    util
