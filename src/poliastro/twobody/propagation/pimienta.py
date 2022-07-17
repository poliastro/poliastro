import sys

from astropy import units as u

from poliastro.core.propagation import pimienta_coe as pimienta_fast
from poliastro.twobody.propagation.enums import PropagatorKind
from poliastro.twobody.states import ClassicalState

from ._compat import OldPropagatorModule

sys.modules[__name__].__class__ = OldPropagatorModule


class PimientaPropagator:
    """Kepler solver for elliptic orbits based on a 15th
    order polynomial with accuracies around 10e-5 for elliptic case and 10e-13
    in the hyperbolic regime.

    Notes
    -----
    This algorithm was developed by Pimienta-Peñalver and John L. Crassidis in
    their paper *Accurate Kepler Equation solver without trascendental function
    evaluations*. Original paper is on Buffalo's UBIR repository: http://hdl.handle.net/10477/50522

    """

    kind = PropagatorKind.ELLIPTIC

    def propagate(self, state, tof):
        state = state.to_classical()

        nu = (
            pimienta_fast(
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
