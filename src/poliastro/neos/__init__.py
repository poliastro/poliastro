"""Code related to NEOs.

Functions related to NEOs and different NASA APIs.
All of them are coded as part of SOCIS 2017 proposal.

Notes
-----

The orbits returned by the functions in this package are in the
:py:class:`~poliastro.frames.HeliocentricEclipticJ2000` frame.

"""

def neos():
    raise NotImplementedError(
        """This module is deprecated from poliastro and a wrapper
        is provided in ~poliastro.twobody.orbit.Orbit for the same.
        For using it, do:

>>> from poliastro.twobody.orbit import Orbit
>>> Orbit.from_sbdb("Neos Name")

instead."""
    )


neows = neos()
