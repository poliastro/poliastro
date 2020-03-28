from enum import Enum, auto

import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astropy.time import Time
from scipy.interpolate import interp1d


class InterpolationMethods(Enum):
    SINC = auto()
    SPLINES = auto()


def build_ephem_interpolant(body, period, t_span, rtol=1e-5):
    """Interpolates ephemerides data

    Parameters
    ----------
    body : Body
        Source body.
    period : ~astropy.units.Quantity
        Orbital period.
    t_span : list(~astropy.units.Quantity)
        Initial and final epochs.
    rtol : float, optional
        Relative tolerance. Controls the number of sampled data points,
        defaults to 1e-5.

    Returns
    -------
    intrp : ~scipy.interpolate.interpolate.interp1d
        Interpolated function.

    """
    h = (period * rtol).to(u.day).value
    t_span = (t_span[0].to(u.day).value, t_span[1].to(u.day).value + 0.01)
    t_values = np.linspace(*t_span, int((t_span[1] - t_span[0]) / h))
    r_values = np.zeros((t_values.shape[0], 3))

    for i, t in enumerate(t_values):
        epoch = Time(t, format="jd", scale="tdb")

        r = get_body_barycentric(body.name, epoch)
        r = (
            ICRS(x=r.x, y=r.y, z=r.z, representation_type=CartesianRepresentation)
            .transform_to(GCRS(obstime=epoch))
            .represent_as(CartesianRepresentation)
        )

        r_values[i] = r.xyz.to(u.km)

    t_values = ((t_values - t_span[0]) * u.day).to(u.s).value
    return interp1d(t_values, r_values, kind="cubic", axis=0, assume_sorted=True)


def _interpolate_splines(epochs, reference_epochs, coordinates, *, kind="cubic"):
    xyz_unit = coordinates.xyz.unit
    d_xyz_unit = coordinates.differentials["s"].d_xyz.unit

    # TODO: Avoid building interpolant every time?
    # Does it have a performance impact?
    interpolant_xyz = interp1d(reference_epochs.jd, coordinates.xyz.value, kind=kind)
    interpolant_d_xyz = interp1d(
        reference_epochs.jd, coordinates.differentials["s"].d_xyz.value, kind=kind,
    )

    result_xyz = interpolant_xyz(epochs.jd) * xyz_unit
    result_d_xyz = interpolant_d_xyz(epochs.jd) * d_xyz_unit

    return CartesianRepresentation(
        result_xyz, differentials=CartesianDifferential(result_d_xyz)
    )


def _interpolate_sinc(epochs, reference_epochs, coordinates):
    def _interp_1d(arr):
        return _sinc_interp(arr, reference_epochs.jd, epochs.jd)

    xyz_unit = coordinates.xyz.unit
    d_xyz_unit = coordinates.differentials["s"].d_xyz.unit

    x = _interp_1d(coordinates.x.value) * xyz_unit
    y = _interp_1d(coordinates.y.value) * xyz_unit
    z = _interp_1d(coordinates.z.value) * xyz_unit

    d_x = _interp_1d(coordinates.differentials["s"].d_x.value) * d_xyz_unit
    d_y = _interp_1d(coordinates.differentials["s"].d_y.value) * d_xyz_unit
    d_z = _interp_1d(coordinates.differentials["s"].d_z.value) * d_xyz_unit

    return CartesianRepresentation(
        x, y, z, differentials=CartesianDifferential(d_x, d_y, d_z)
    )


# Taken from https://gist.github.com/endolith/1297227
def _sinc_interp(x, s, u):
    """Interpolates x, sampled at "s" instants, at "u" instants.

    """
    if len(x) != len(s):
        raise ValueError("x and s must be the same length")

    # Find the period and assume it's constant
    T = s[1] - s[0]

    sincM = np.tile(u, (len(s), 1)) - np.tile(s[:, np.newaxis], (1, len(u)))
    y = np.dot(x, np.sinc(sincM / T))

    return y


_INTERPOLATION_MAPPING = {
    InterpolationMethods.SINC: _interpolate_sinc,
    InterpolationMethods.SPLINES: _interpolate_splines,
}


class Ephem:
    def __init__(self, coordinates, epochs):
        self._epochs = epochs
        self._coordinates = coordinates

    @property
    def epochs(self):
        return self._epochs

    @classmethod
    def from_body(cls, body, epochs):
        """Return `Ephem` from a body ephemerides at certain epochs.

        Parameters
        ----------
        body: ~poliastro.bodies.Body
            Body.
        epochs: ~astropy.time.Time
            Epochs to sample the body positions.

        """
        r, v = get_body_barycentric_posvel(body.name, epochs)
        coordinates = r.with_differentials(v.represent_as(CartesianDifferential))

        return cls(coordinates, epochs)

    def sample(self, epochs=None, *, method=InterpolationMethods.SPLINES, **kwargs):
        if epochs is None:
            return self._coordinates

        # TODO: Proper type annotation
        coordinates = _INTERPOLATION_MAPPING[method](
            epochs, self.epochs, self._coordinates, **kwargs
        )  # type: ignore

        return coordinates
