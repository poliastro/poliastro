"""Coordinate frames definitions.

"""
import numpy as np

from astropy import _erfa
from astropy import units as u
from astropy.coordinates import (
    get_body_barycentric, frame_transform_graph,
    BaseEclipticFrame, BaseRADecFrame,
    ICRS, HCRS, GCRS,
    TimeAttribute,
    AffineTransform, FunctionTransformWithFiniteDifference,
    UnitSphericalRepresentation,
)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME, get_jd12
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_transpose

from poliastro.constants import J2000
from poliastro.bodies import Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto


class HeliocentricEclipticJ2000(BaseEclipticFrame):
    """
    Heliocentric ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the mean equinox of J2000 and the xy-plane in the plane of the
    ecliptic of J2000 (according to the IAU 1976/1980 obliquity model).

    """
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


def _ecliptic_rotation_matrix():
    jd1, jd2 = get_jd12(J2000, J2000.scale)
    obl = _erfa.obl80(jd1, jd2) * u.radian
    assert obl.to(u.arcsec).value == 84381.448
    return rotation_matrix(obl, 'x')


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 ICRS, HeliocentricEclipticJ2000,
                                 finite_difference_frameattr_name='obstime')
def _icrs_to_ecliptic(from_coo, to_frame):
    # get barycentric sun coordinate
    bary_sun_pos = get_body_barycentric('sun', to_frame.obstime)

    # offset to heliocentric
    heliocart = from_coo.cartesian - bary_sun_pos

    # now compute the matrix to precess to the right orientation
    rmat = _ecliptic_rotation_matrix()

    newrepr = heliocart.transform(rmat)
    return to_frame.realize_frame(newrepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference,
                                 HeliocentricEclipticJ2000, ICRS,
                                 finite_difference_frameattr_name='obstime')
def _ecliptic_to_icrs(from_coo, to_frame):
    # first un-precess from ecliptic to ICRS orientation
    rmat = _ecliptic_rotation_matrix()
    intermed_repr = from_coo.cartesian.transform(matrix_transpose(rmat))

    # now offset back to barycentric, which is the correct center for ICRS
    # get barycentric sun coordinate
    bary_sun_pos = get_body_barycentric('sun', from_coo.obstime)

    newrepr = intermed_repr + bary_sun_pos
    return to_frame.realize_frame(newrepr)


# ---
# Planetary frames parallel to ICRS
# Taken from Astropy HCRS
# ---

_NEED_ORIGIN_HINT = ("The input {0} coordinates do not have length units. This "
                     "probably means you created coordinates with lat/lon but "
                     "no distance.  PlanetaryICRS<->ICRS transforms cannot "
                     "function in this case because there is an origin shift.")


class _PlanetaryICRS(BaseRADecFrame):
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(AffineTransform, cls, ICRS)(cls.to_icrs)
        frame_transform_graph.transform(AffineTransform, ICRS, cls)(cls.from_icrs)
        frame_transform_graph.transform(FunctionTransformWithFiniteDifference, cls, cls)(cls.self_transform)

        return super().__new__(cls)

    @staticmethod
    def to_icrs(planet_coo, _):
        # this is just an origin translation so without a distance it cannot go ahead
        if isinstance(planet_coo.data, UnitSphericalRepresentation):
            raise u.UnitsError(_NEED_ORIGIN_HINT.format(planet_coo.__class__.__name__))

        if planet_coo.data.differentials:
            from astropy.coordinates.solar_system import get_body_barycentric_posvel
            bary_sun_pos, bary_sun_vel = get_body_barycentric_posvel(planet_coo.body.name,
                                                                     planet_coo.obstime)
            bary_sun_pos = bary_sun_pos.with_differentials(bary_sun_vel)

        else:
            from astropy.coordinates.solar_system import get_body_barycentric
            bary_sun_pos = get_body_barycentric(planet_coo.body.name, planet_coo.obstime)
            bary_sun_vel = None

        return None, bary_sun_pos

    @staticmethod
    def from_icrs(icrs_coo, planet_frame):
        # this is just an origin translation so without a distance it cannot go ahead
        if isinstance(icrs_coo.data, UnitSphericalRepresentation):
            raise u.UnitsError(_NEED_ORIGIN_HINT.format(icrs_coo.__class__.__name__))

        if icrs_coo.data.differentials:
            from astropy.coordinates.solar_system import get_body_barycentric_posvel
            bary_sun_pos, bary_sun_vel = get_body_barycentric_posvel(planet_frame.body.name,
                                                                     planet_frame.obstime)
            bary_sun_pos = -bary_sun_pos.with_differentials(-bary_sun_vel)

        else:
            from astropy.coordinates.solar_system import get_body_barycentric
            bary_sun_pos = -get_body_barycentric(planet_frame.body.name, planet_frame.obstime)
            bary_sun_vel = None

        return None, bary_sun_pos

    @staticmethod
    def self_transform(from_coo, to_frame):
        if np.all(from_coo.obstime == to_frame.obstime):
            return to_frame.realize_frame(from_coo.data)
        else:
            # like CIRS, we do this self-transform via ICRS
            return from_coo.transform_to(ICRS).transform_to(to_frame)


class MercuryICRS(_PlanetaryICRS):
    body = Mercury


class VenusICRS(_PlanetaryICRS):
    body = Venus


class MarsICRS(_PlanetaryICRS):
    body = Mars


class JupiterICRS(_PlanetaryICRS):
    body = Jupiter


class SaturnICRS(_PlanetaryICRS):
    body = Saturn


class UranusICRS(_PlanetaryICRS):
    body = Uranus


class NeptuneICRS(_PlanetaryICRS):
    body = Neptune


class PlutoICRS(_PlanetaryICRS):
    body = Pluto


INERTIAL_FRAME_MAPPING = {
    Sun: (HCRS,),
    Mercury: (MercuryICRS,),
    Venus: (VenusICRS,),
    Earth: (GCRS,),
    Mars: (MarsICRS,),
    Jupiter: (JupiterICRS,),
    Saturn: (SaturnICRS,),
    Uranus: (UranusICRS,),
    Neptune: (NeptuneICRS,),
    Pluto: (PlutoICRS,)
}
