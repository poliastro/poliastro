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

from .bodies import Earth, Jupiter, Mars, Mercury, Neptune, Pluto, Saturn, Uranus, Venus
from .constants import J2000
from .frames import Planes
from .twobody.angles import M_to_nu
from .twobody.orbit import Orbit

# Source: https://ssd.jpl.nasa.gov/?planet_pos
MEAN_ELEMENTS_1800AD_2050AD = {
    Mercury: {
        "a": (0.38709927, 0.00000037),
        "ecc": (0.20563593, 0.00001906),
        "inc": (7.00497902, -0.00594749),
        "meanlon": (252.25032350, 149472.67411175),
        "lonper": (77.45779628, 0.16047689),
        "node": (48.33076593, -0.12534081),
    },
    Venus: {
        "a": (0.72333566, 3.9e-06),
        "ecc": (0.00677672, -4.107e-05),
        "inc": (3.39467605, -0.0007889),
        "meanlon": (181.9790995, 58517.81538729),
        "lonper": (131.60246718, 0.00268329),
        "node": (76.67984255, -0.27769418),
    },
    Earth: {
        "a": (1.00000261, 5.62e-06),
        "ecc": (0.01671123, -4.392e-05),
        "inc": (-1.531e-05, -0.01294668),
        "meanlon": (100.46457166, 35999.37244981),
        "lonper": (102.93768193, 0.32327364),
        "node": (0.0, 0.0),
    },
    Mars: {
        "a": (1.52371034, 1.847e-05),
        "ecc": (0.0933941, 7.882e-05),
        "inc": (1.84969142, -0.00813131),
        "meanlon": (-4.55343205, 19140.30268499),
        "lonper": (-23.94362959, 0.44441088),
        "node": (49.55953891, -0.29257343),
    },
    Jupiter: {
        "a": (5.202887, -0.00011607),
        "ecc": (0.04838624, -0.00013253),
        "inc": (1.30439695, -0.00183714),
        "meanlon": (34.39644051, 3034.74612775),
        "lonper": (14.72847983, 0.21252668),
        "node": (100.47390909, 0.20469106),
    },
    Saturn: {
        "a": (9.53667594, -0.0012506),
        "ecc": (0.05386179, -0.00050991),
        "inc": (2.48599187, 0.00193609),
        "meanlon": (49.95424423, 1222.49362201),
        "lonper": (92.59887831, -0.41897216),
        "node": (113.66242448, -0.28867794),
    },
    Uranus: {
        "a": (19.18916464, -0.00196176),
        "ecc": (0.04725744, -4.397e-05),
        "inc": (0.77263783, -0.00242939),
        "meanlon": (313.23810451, 428.48202785),
        "lonper": (170.9542763, 0.40805281),
        "node": (74.01692503, 0.04240589),
    },
    Neptune: {
        "a": (30.06992276, 0.00026291),
        "ecc": (0.00859048, 5.105e-05),
        "inc": (1.77004347, 0.00035372),
        "meanlon": (-55.12002969, 218.45945325),
        "lonper": (44.96476227, -0.32241464),
        "node": (131.78422574, -0.00508664),
    },
    Pluto: {
        "a": (39.48211675, -0.00031596),
        "ecc": (0.2488273, 5.17e-05),
        "inc": (17.14001206, 4.818e-05),
        "meanlon": (238.92903833, 145.20780515),
        "lonper": (224.06891629, -0.04062942),
        "node": (110.30393684, -0.01183482),
    },
}


def _get_element(element_0, element_dot, epoch=J2000):
    T_cent = (epoch.jd - J2000.jd) / 36525
    return element_0 + element_dot * T_cent


def _get_elements(body, epoch=J2000):
    elements = MEAN_ELEMENTS_1800AD_2050AD[body]
    return {
        name: _get_element(elem[0], elem[1], epoch) for name, elem in elements.items()
    }


def get_mean_orbit(body, epoch=J2000):
    elements = _get_elements(body, epoch)

    a = elements["a"]
    e = elements["ecc"]
    inc = elements["inc"]
    node = elements["node"]

    if inc < 0:
        inc = -inc
        node = (node + 180) % 360

    argp = elements["lonper"] - node
    M = elements["meanlon"] - elements["lonper"]

    return Orbit.from_classical(
        body.parent,
        a * u.au,
        e * u.one,
        inc * u.deg,
        node * u.deg,
        argp * u.deg,
        M_to_nu(M * u.deg, e * u.one),
        epoch,
        plane=Planes.EARTH_ECLIPTIC,
    )


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
