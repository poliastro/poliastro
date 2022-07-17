import sys

import numpy as np
from astropy import units as u

from poliastro.core.propagation.farnocchia import (
    farnocchia_coe as farnocchia_coe_fast,
    farnocchia_rv as farnocchia_rv_fast,
)
from poliastro.twobody.propagation.enums import PropagatorKind
from poliastro.twobody.states import ClassicalState

from ._compat import OldPropagatorModule

sys.modules[__name__].__class__ = OldPropagatorModule


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
            farnocchia_coe_fast(
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

    def propagate_many(self, state, tofs):
        state = state.to_vectors()
        k = state.attractor.k.to_value(u.km**3 / u.s**2)
        rv0 = state.to_value()

        # TODO: This should probably return a ClassicalStateArray instead,
        # see discussion at https://github.com/poliastro/poliastro/pull/1492
        results = np.array(
            [farnocchia_rv_fast(k, *rv0, tof) for tof in tofs.to_value(u.s)]
        )
        return (
            results[:, 0] << u.km,
            results[:, 1] << (u.km / u.s),
        )
