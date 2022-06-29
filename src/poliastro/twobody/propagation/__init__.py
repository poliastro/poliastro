"""The following script holds the different high level functions for the
different propagators available at poliastro:

+-------------+------------+-----------------+-----------------+
|  Propagator | Elliptical |    Parabolic    |    Hyperbolic   |
+-------------+------------+-----------------+-----------------+
|  farnocchia |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|   vallado   |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|   mikkola   |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|   markley   |      ✓     |        x        |        x        |
+-------------+------------+-----------------+-----------------+
|   pimienta  |      ✓     |        ✓        |        x        |
+-------------+------------+-----------------+-----------------+
|   gooding   |      ✓     |        x        |        x        |
+-------------+------------+-----------------+-----------------+
|    danby    |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|    cowell   |      ✓     |        ✓        |        ✓        |
+-------------+------------+-----------------+-----------------+
|  recseries  |      ✓     |        x        |        x        |
+-------------+------------+-----------------+-----------------+

"""
# FIXME: Make submodules private, to avoid "module not callable" errors?
from poliastro.twobody.propagation.cowell import CowellPropagator
from poliastro.twobody.propagation.danby import DanbyPropagator
from poliastro.twobody.propagation.enums import PropagatorKind
from poliastro.twobody.propagation.farnocchia import FarnocchiaPropagator
from poliastro.twobody.propagation.gooding import GoodingPropagator
from poliastro.twobody.propagation.markley import MarkleyPropagator
from poliastro.twobody.propagation.mikkola import MikkolaPropagator
from poliastro.twobody.propagation.pimienta import PimientaPropagator
from poliastro.twobody.propagation.recseries import RecseriesPropagator
from poliastro.twobody.propagation.vallado import ValladoPropagator


# FIXME: Make `propagate` private API?
def propagate(
    state,
    time_of_flight,
    *,
    method=FarnocchiaPropagator(),
):
    """Propagate an orbit some time and return the result.

    Parameters
    ----------
    state : ~poliastro.twobody.states.BaseState
        State object to propagate.
    time_of_flight : ~astropy.time.TimeDelta
        Time of propagation.
    method : callable, optional
        Propagation method, default to farnocchia.

    Returns
    -------
    astropy.coordinates.CartesianRepresentation
        Propagation coordinates.
    """

    # Check if propagator fulfills orbit requirements
    # Note there's a potential conversion here purely for convenience that could be skipped
    ecc = state.to_classical().ecc
    if ecc < 1.0 and not (method.kind & PropagatorKind.ELLIPTIC):
        raise ValueError(
            "Can not use an parabolic/hyperbolic propagator for elliptical/circular orbits."
        )
    elif ecc == 1.0 and not (method.kind & PropagatorKind.PARABOLIC):
        raise ValueError(
            "Can not use an elliptic/hyperbolic propagator for parabolic orbits."
        )
    elif ecc > 1.0 and not (method.kind & PropagatorKind.HYPERBOLIC):
        raise ValueError(
            "Can not use an elliptic/parabolic propagator for hyperbolic orbits."
        )

    if time_of_flight.ndim != 0:
        raise ValueError(
            "propagate only accepts scalar values for time_of_flight"
        )

    new_state = method.propagate(state, time_of_flight)

    return new_state


ALL_PROPAGATORS = [
    CowellPropagator,
    DanbyPropagator,
    FarnocchiaPropagator,
    GoodingPropagator,
    MarkleyPropagator,
    MikkolaPropagator,
    PimientaPropagator,
    RecseriesPropagator,
    ValladoPropagator,
]
ELLIPTIC_PROPAGATORS = [
    propagator
    for propagator in ALL_PROPAGATORS
    if propagator.kind & PropagatorKind.ELLIPTIC
]
PARABOLIC_PROPAGATORS = [
    propagator
    for propagator in ALL_PROPAGATORS
    if propagator.kind & PropagatorKind.PARABOLIC
]
HYPERBOLIC_PROPAGATORS = [
    propagator
    for propagator in ALL_PROPAGATORS
    if propagator.kind & PropagatorKind.HYPERBOLIC
]
