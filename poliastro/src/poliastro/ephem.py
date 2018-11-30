import numpy as np
from scipy.interpolate import interp1d

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, ICRS, GCRS, CartesianRepresentation


def build_ephem_interpolant(body, period, t_span, rtol=1e-5):
    h = (period * rtol).to(u.day).value
    t_span = ((t_span[0].to(u.day).value, t_span[1].to(u.day).value + 0.01))
    t_values = np.linspace(*t_span, int((t_span[1] - t_span[0]) / h))
    r_values = np.zeros((t_values.shape[0], 3))

    for i, t in enumerate(t_values):
        epoch = Time(t, format='jd', scale='tdb')

        r = get_body_barycentric(body.name, epoch)
        r = (ICRS(x=r.x, y=r.y, z=r.z, representation_type=CartesianRepresentation)
             .transform_to(GCRS(obstime=epoch))
             .represent_as(CartesianRepresentation)
             )

        r_values[i] = r.xyz.to(u.km)

    t_values = ((t_values - t_span[0]) * u.day).to(u.s).value
    return interp1d(t_values, r_values, kind='cubic', axis=0, assume_sorted=True)
