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
