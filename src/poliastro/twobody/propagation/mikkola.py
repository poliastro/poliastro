import sys

from astropy import units as u

from poliastro.core.propagation import mikkola_coe as mikkola_fast
from poliastro.twobody.propagation.enums import PropagatorKind
from poliastro.twobody.states import ClassicalState

from ._compat import OldPropagatorModule

sys.modules[__name__].__class__ = OldPropagatorModule


class MikkolaPropagator:
    """Solves Kepler Equation by a cubic approximation.

    Notes
    -----
    This method was derived by Seppo Mikola in his paper *A Cubic Approximation
    For Kepler's Equation* with DOI: https://doi.org/10.1007/BF01235850

    """

    kind = PropagatorKind.ELLIPTIC | PropagatorKind.HYPERBOLIC

    def propagate(self, state, tof):
        state = state.to_classical()

        nu = (
            mikkola_fast(
                state.attractor.k.to_value(u.km**3 / u.s**2),
                *state.to_value(),
                tof.to_value(u.s)
            )
            << u.rad
        )

        new_state = ClassicalState(
            state.attractor, state.to_tuple()[:5] + (nu,), state.plane
        )
        return new_state
