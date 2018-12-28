Neos module
===========

NEO stans for near-Earth object. These bodies are generally asteroids or coments such
that because of their orbits could be close into Earth's surroundings. NEOs can be classified
into:

    * Meteroids: with a diameter lower than 50 meters
    * Near Earth Comets (NEC): with orbital period lower than 200 years
    * Near Earth Asteroid (NEA): most of NEO objects are in this category

:py:mod:`poliastro.neos` was developed as a result for the SOCIS 2017 Edition. It works by
sending a request to  and `NASANeo`_ webpage downloading all the parameters that define the orbit of the body.
This means that Internet conection is required to accomplish this task.

.. _NASANeo: https://api.nasa.gov/api.html#NeoWS


.. toctree::
    :maxdepth: 1

    dastcom5
    dastcom5_parameters
    neows
