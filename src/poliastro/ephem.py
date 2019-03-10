import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    CartesianRepresentation,
    get_body_barycentric,
)
from astropy.time import Time
from scipy.interpolate import interp1d


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
