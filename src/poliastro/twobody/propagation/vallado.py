import sys

import numpy as np
from astropy import units as u

from poliastro.core.propagation import vallado as vallado_fast
from poliastro.twobody.propagation.enums import PropagatorKind
from poliastro.twobody.states import RVState

from ._compat import OldPropagatorModule

sys.modules[__name__].__class__ = OldPropagatorModule


def vallado(k, r0, v0, tof, *, numiter):
    # Compute Lagrange coefficients
    f, g, fdot, gdot = vallado_fast(k, r0, v0, tof, numiter)

    assert (
        np.abs(f * gdot - fdot * g - 1) < 1e-5
    ), "Internal error, solution is not consistent"  # Fixed tolerance

    # Return position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v


class ValladoPropagator:
    """Propagates Keplerian orbit using Vallado's method.

    Notes
    -----
    This algorithm is based on Vallado implementation, and does basic Newton
    iteration on the Kepler equation written using universal variables. Battin
    claims his algorithm uses the same amount of memory but is between 40 %
    and 85 % faster.

    """

    kind = (
        PropagatorKind.ELLIPTIC
        | PropagatorKind.PARABOLIC
        | PropagatorKind.HYPERBOLIC
    )

    def __init__(self, numiter=350):
        self._numiter = numiter

    def propagate(self, state, tof):
        state = state.to_vectors()

        r_raw, v_raw = vallado(
            state.attractor.k.to_value(u.km**3 / u.s**2),
            *state.to_value(),
            tof.to_value(u.s),
            numiter=self._numiter,
        )
        r = r_raw << u.km
        v = v_raw << (u.km / u.s)

        new_state = RVState(state.attractor, (r, v), state.plane)
        return new_state
