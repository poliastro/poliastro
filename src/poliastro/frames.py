"""Coordinate frames definitions.

"""
from astropy import _erfa
from astropy import units as u
from astropy.coordinates import (
    get_body_barycentric, frame_transform_graph,
    BaseEclipticFrame, ICRS,
    TimeAttribute, FunctionTransformWithFiniteDifference,
)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME, get_jd12
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_transpose

from poliastro.constants import J2000


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

