# coding: utf-8
import numpy as np

from astropy import units as u

from astropy import time

from poliastro.twobody.propagation import kepler
from poliastro.twobody.core import StateFactory

J2000 = time.Time("J2000", scale='utc')


def propagate(orbit, time_of_flight, rtol=1e-10):
    """Propagate this `State` some `time` and return the result.

    """
    r, v = kepler(orbit.state.attractor.k.to(u.km ** 3 / u.s ** 2).value,
                  orbit.state.r.to(u.km).value, orbit.state.v.to(u.km / u.s).value,
                  time_of_flight.to(u.s).value,
                  rtol=rtol)
    new_ss = StateFactory.from_vectors(orbit.state.attractor, r * u.km, v * u.km / u.s)
    return Orbit(new_ss, orbit.epoch + time_of_flight)


class Orbit(object):
    def __init__(self, state, epoch=J2000):
        self.state = state
        self.epoch = epoch

    def propagate(self, time_of_flight, rtol=1e-10):
        """Propagate this `State` some `time` and return the result.

        """
        return propagate(self, time_of_flight, rtol)

    def apply_maneuver(self, maneuver, intermediate=False):
        """Returns resulting State after applying maneuver to self.

        Optionally return intermediate states (default to False).

        Parameters
        ----------
        maneuver : Maneuver
            Maneuver to apply.
        intermediate : bool, optional
            Return intermediate states, default to False.

        """
        orbit_new = self  # Initialize
        states = []
        attractor = self.state.attractor
        for delta_t, delta_v in maneuver:
            if not delta_t == 0 * u.s:
                orbit_new = self.propagate(time_of_flight=delta_t)
            r, v = orbit_new.state.rv()
            vnew = v + delta_v
            ss_new = StateFactory.from_vectors(attractor, r, vnew)
            orbit_new = Orbit(ss_new, orbit_new.epoch)
            states.append(orbit_new)
        if intermediate:
            res = states
        else:
            res = orbit_new
        return res
