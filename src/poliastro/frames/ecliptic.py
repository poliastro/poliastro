import erfa
import numpy as np
from astropy import units as u
from astropy.coordinates import (
    BaseEclipticFrame,
    CartesianRepresentation,
    DynamicMatrixTransform,
    GeocentricMeanEcliptic as _GeocentricMeanEcliptic,
    HeliocentricEclipticIAU76 as _HeliocentricEclipticJ2000,
    TimeAttribute,
    UnitSphericalRepresentation,
    frame_transform_graph,
    get_body,
)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME, get_jd12
from astropy.coordinates.matrix_utilities import (
    matrix_product,
    matrix_transpose,
    rotation_matrix,
)

from .equatorial import GCRS

__all__ = [
    "GeocentricSolarEcliptic",
    "GeocentricMeanEcliptic",
    "HeliocentricEclipticJ2000",
]

# HACK: sphinx-autoapi variable definition
GeocentricMeanEcliptic = _GeocentricMeanEcliptic
HeliocentricEclipticJ2000 = _HeliocentricEclipticJ2000


class GeocentricSolarEcliptic(BaseEclipticFrame):
    """
    This system has its X axis towards the Sun and its Z axis perpendicular to
    the plane of the Earth's orbit around the Sun (positive North). This system
    is fixed with respect to the Earth-Sun line. It is convenient for specifying
    magnetospheric boundaries. It has also been widely adopted as the system for
    representing vector quantities in space physics databases.

    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


@frame_transform_graph.transform(DynamicMatrixTransform, GCRS, GeocentricSolarEcliptic)
def gcrs_to_geosolarecliptic(gcrs_coo, to_frame):

    if not to_frame.obstime.isscalar:
        raise ValueError(
            "To perform this transformation the obstime Attribute must be a scalar."
        )

    _earth_orbit_perpen_point_gcrs = UnitSphericalRepresentation(
        lon=0 * u.deg, lat=(90 * u.deg - _obliquity_rotation_value(to_frame.obstime))
    )

    _earth_detilt_matrix = _make_rotation_matrix_from_reprs(
        _earth_orbit_perpen_point_gcrs, CartesianRepresentation(0, 0, 1)
    )

    sun_pos_gcrs = get_body("sun", to_frame.obstime).cartesian
    earth_pos_gcrs = get_body("earth", to_frame.obstime).cartesian
    sun_earth = sun_pos_gcrs - earth_pos_gcrs

    sun_earth_detilt = sun_earth.transform(_earth_detilt_matrix)

    # Earth-Sun Line in Geocentric Solar Ecliptic Frame
    x_axis = CartesianRepresentation(1, 0, 0)

    rot_matrix = _make_rotation_matrix_from_reprs(sun_earth_detilt, x_axis)

    return matrix_product(rot_matrix, _earth_detilt_matrix)


@frame_transform_graph.transform(DynamicMatrixTransform, GeocentricSolarEcliptic, GCRS)
def geosolarecliptic_to_gcrs(from_coo, gcrs_frame):
    return matrix_transpose(gcrs_to_geosolarecliptic(gcrs_frame, from_coo))


def _obliquity_rotation_value(equinox):
    """
    Function to calculate obliquity of the earth.
    This uses obl06 of erfa.
    """
    jd1, jd2 = get_jd12(equinox, "tt")
    obl = erfa.obl06(jd1, jd2) * u.radian
    return obl.to(u.deg)


def _make_rotation_matrix_from_reprs(start_representation, end_representation):
    """
    Return the matrix for the direct rotation from one representation to a second representation.
    The representations need not be normalized first.
    """
    A = start_representation.to_cartesian()
    B = end_representation.to_cartesian()
    rotation_axis = A.cross(B)
    rotation_angle = -np.arccos(
        A.dot(B) / (A.norm() * B.norm())
    )  # Negation is required

    # This line works around some input/output quirks of Astropy's rotation_matrix()
    matrix = np.array(rotation_matrix(rotation_angle, rotation_axis.xyz.value.tolist()))
    return matrix
