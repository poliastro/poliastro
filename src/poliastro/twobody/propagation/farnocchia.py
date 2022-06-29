from astropy import units as u

from poliastro.core.propagation.farnocchia import (
    farnocchia_coe as farnocchia_fast,
)
from poliastro.twobody.propagation.enums import PropagatorKind
from poliastro.twobody.states import ClassicalState


class FarnocchiaPropagator:
    r"""Propagates orbit using Farnocchia's method.

    Notes
    -----
    This method takes initial :math:`\vec{r}, \vec{v}`, calculates classical orbit parameters,
    increases mean anomaly and performs inverse transformation to get final :math:`\vec{r}, \vec{v}`
    The logic is based on formulae (4), (6) and (7) from http://dx.doi.org/10.1007/s10569-013-9476-9

    """

    kind = (
        PropagatorKind.ELLIPTIC
        | PropagatorKind.PARABOLIC
        | PropagatorKind.HYPERBOLIC
    )

    def propagate(self, state, tof):
        state = state.to_classical()

        nu = (
            farnocchia_fast(
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
