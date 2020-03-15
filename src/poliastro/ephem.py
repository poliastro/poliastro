import numpy as np
from astropy import _erfa as erfa, units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    CartesianRepresentation,
    get_body_barycentric,
)
from astropy.coordinates.solar_system import PLAN94_BODY_NAME_TO_PLANET_INDEX
from astropy.time import Time
from scipy.interpolate import interp1d

from .constants import J2000
from .frames import Planes
from .twobody.states import RVState


def get_mean_elements(body, epoch=J2000):
    """Get mean elements of body.

    """
    # Internally, erfa.plan94 has a table of classical elements,
    # which are then converted to position and velocity...
    # So we are reverting the conversion in the last line
    # This way at least we avoid copy pasting the values
    name = body.name.lower()
    if name in ("earth", "moon"):
        name = "earth-moon-barycenter"

    body_index = PLAN94_BODY_NAME_TO_PLANET_INDEX[name]
    body_pv_helio = erfa.plan94(epoch.jd1, epoch.jd2, body_index)

    return RVState(
        body.parent,
        body_pv_helio["p"] * u.au,
        body_pv_helio["v"] * u.au / u.day,
        plane=Planes.EARTH_ECLIPTIC,
    ).to_classical()


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
