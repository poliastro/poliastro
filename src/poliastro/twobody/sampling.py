import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation

from poliastro.twobody.angles import E_to_nu, nu_to_E
from poliastro.twobody.elements import coe2rv_many, hyp_nu_limit, t_p
from poliastro.twobody.propagation import FarnocchiaPropagator
from poliastro.util import alinspace, wrap_angle


@u.quantity_input(ecc=u.one, min_nu=u.rad, max_nu=u.rad)
def sample_closed(ecc, min_nu, max_nu=None, num_values=100):
    """Sample a closed orbit.

    If ``max_nu`` is given, the sampling interval will go
    from the minimum to the maximum true anomaly in the direction of the orbit.
    If not given, it will do a full revolution starting in the minimum true anomaly.

    Notes
    -----
    First sample the eccentric anomaly uniformly,
    then transform into true anomaly
    to minimize error in the apocenter,
    see https://apps.dtic.mil/dtic/tr/fulltext/u2/a605040.pdf

    """
    # Because how nu_to_E works, we don't need to wrap the angle here!
    # It will do the right thing
    min_E = nu_to_E(min_nu, ecc)

    # This linspace will always increase positively,
    # even though it might contain out of range values
    E_values = alinspace(
        min_E,
        nu_to_E(max_nu, ecc) if max_nu is not None else None,
        num=num_values,
    )

    # Because how E_to_nu works, we don't need to wrap the angles here!
    # It will do the right thing
    nu_values = E_to_nu(E_values, ecc)

    # We wrap the angles on return
    return (nu_values + np.pi * u.rad) % (2 * np.pi * u.rad) - np.pi * u.rad


@u.quantity_input(min_nu=u.rad, ecc=u.one, max_nu=u.rad, nu_limit=u.rad)
def sample_open(
    ecc, min_nu=None, max_nu=None, num_values=100, *, nu_limit=None
):
    """Sample an open orbit.

    Notes
    -----
    Uniform sampling on true anomaly in the absence of a better method.
    Minimum and maximum anomaly must be within limits,
    which are computed from the eccentricity if not given.

    """
    if nu_limit is None:
        nu_limit = hyp_nu_limit(ecc)

    min_nu = min_nu if min_nu is not None else -nu_limit
    max_nu = max_nu if max_nu is not None else nu_limit

    if not (-nu_limit) <= min_nu < max_nu <= nu_limit:
        raise ValueError("Anomaly values out of range")

    nu_values = alinspace(min_nu, max_nu, num=num_values)

    # No need to wrap the angles
    return nu_values


class SamplingStrategy:
    def sample(self, orbit):
        raise NotImplementedError


class EpochsArray(SamplingStrategy):
    def __init__(self, epochs, method=FarnocchiaPropagator()):
        self._epochs = epochs
        self._method = method

    def sample(self, orbit):
        times_of_flight = self._epochs - orbit.epoch
        # TODO: Make state public?
        rr, vv = self._method.propagate_many(orbit._state, times_of_flight)

        # TODO: For now we're returning CartesianRepresentation
        # because that's what Ephem objects expect,
        # but we could probably use Classical/RVStateArray instead when available.
        # However, we are also returning the epochs
        # (since computing them here is more efficient than doing it from the outside)
        # but there are open questions around StateArrays and epochs.
        # See discussion at https://github.com/poliastro/poliastro/pull/1492
        cartesian = CartesianRepresentation(
            rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
        )
        return cartesian, self._epochs


class TrueAnomalyBounds(SamplingStrategy):
    def __init__(
        self, min_nu=None, max_nu=None, num_values=100, hyp_r_factor=3.0
    ):
        self._min_nu = min_nu
        self._max_nu = max_nu
        self._hyp_r_factor = hyp_r_factor
        self._num_values = num_values

    def sample(self, orbit):
        if orbit.ecc < 1:
            nu_values = sample_closed(
                orbit.ecc,
                self._min_nu if self._min_nu is not None else orbit.nu,
                self._max_nu,
                self._num_values,
            )
        else:
            # Select a sensible limiting value for non-closed orbits
            # This corresponds to max(r = 3p, r = self.r)
            # We have to wrap nu in [-180, 180) to compare it with the output of
            # the arc cosine, which is in the range [0, 180)
            # Start from -nu_limit
            nu_limit = max(
                hyp_nu_limit(orbit.ecc, self._hyp_r_factor),
                abs(wrap_angle(orbit.nu)),
            )

            # Perform actual sampling
            nu_values = sample_open(
                orbit.ecc,
                self._min_nu,
                self._max_nu,
                self._num_values,
                nu_limit=nu_limit,
            )

        # Weird units roundtrip because list of quantities != array quantity
        delta_ts = [
            t_p(nu, orbit.ecc, orbit.attractor.k, orbit.r_p).to_value(u.s)
            for nu in nu_values
        ]
        # Unwrap time increments to return monotonic increasing epochs
        # Notice that astropy.units does not support the period kwarg for np.unwrap
        delta_ts = (
            np.unwrap(delta_ts, period=orbit.period.to_value(u.s)) << u.s
        )
        epochs = orbit.epoch + (orbit.t_p + delta_ts)

        n = nu_values.shape[0]
        rr, vv = coe2rv_many(
            np.tile(orbit.attractor.k, n),
            np.tile(orbit.p, n),
            np.tile(orbit.ecc, n),
            np.tile(orbit.inc, n),
            np.tile(orbit.raan, n),
            np.tile(orbit.argp, n),
            nu_values,
        )

        # TODO: For now we're returning CartesianRepresentation
        # because that's what Ephem objects expect,
        # but we could probably use Classical/RVStateArray instead when available.
        # However, we are also returning the epochs
        # (since computing them here is more efficient than doing it from the outside)
        # but there are open questions around StateArrays and epochs.
        # See discussion at https://github.com/poliastro/poliastro/pull/1492
        cartesian = CartesianRepresentation(
            rr, differentials=CartesianDifferential(vv, xyz_axis=1), xyz_axis=1
        )
        return cartesian, epochs


class EpochBounds(SamplingStrategy):
    def __init__(self, min_epoch=None, max_epoch=None, num_values=100):
        self._min_epoch = min_epoch
        self._max_epoch = max_epoch
        self._num_values = num_values

    def sample(self, orbit):
        if self._min_epoch is None:
            min_nu = orbit.nu
        else:
            min_nu = orbit.propagate(self._min_epoch).nu

        if self._max_epoch is None:
            max_nu = None
        else:
            max_nu = orbit.propagate(self._max_epoch).nu

        return TrueAnomalyBounds(
            min_nu=min_nu, max_nu=max_nu, num_values=self._num_values
        ).sample(orbit)
