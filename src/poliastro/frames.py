"""Coordinate frames definitions.

"""
from enum import Enum
from typing import Dict

import numpy as np
from astropy import _erfa as erfa, units as u
from astropy.coordinates import (
    GCRS,
    ICRS,
    AffineTransform,
    BaseEclipticFrame,
    BaseRADecFrame,
    CartesianDifferential,
    CartesianRepresentation,
    FunctionTransformWithFiniteDifference,
    HeliocentricEclipticIAU76 as HeliocentricEclipticJ2000,
    TimeAttribute,
    UnitSphericalRepresentation,
    frame_transform_graph,
    get_body,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astropy.coordinates.baseframe import FrameMeta
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME, get_jd12
from astropy.coordinates.matrix_utilities import (
    matrix_product,
    matrix_transpose,
    rotation_matrix,
)
from astropy.coordinates.transformations import DynamicMatrixTransform

from poliastro.bodies import (
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
    _Body,
)
from poliastro.constants import J2000


class Planes(Enum):
    EARTH_EQUATOR = "Earth mean Equator and Equinox of epoch (J2000.0)"
    EARTH_ECLIPTIC = "Earth mean Ecliptic and Equinox of epoch (J2000.0)"
    # BODY_EQUATOR = 'Body mean Equator and node of date'  # TODO: Implement proper conversions


# ---
# Planetary frames parallel to ICRS
# Taken from Astropy HCRS
# ---

_NEED_ORIGIN_HINT = (
    "The input {0} coordinates do not have length units. This "
    "probably means you created coordinates with lat/lon but "
    "no distance.  PlanetaryICRS<->ICRS transforms cannot "
    "function in this case because there is an origin shift."
)


class _PlanetaryICRS(BaseRADecFrame):
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(AffineTransform, cls, ICRS)(cls.to_icrs)
        frame_transform_graph.transform(AffineTransform, ICRS, cls)(cls.from_icrs)
        frame_transform_graph.transform(
            FunctionTransformWithFiniteDifference, cls, cls
        )(cls.self_transform)

        return super().__new__(cls)

    @staticmethod
    def to_icrs(planet_coo, _):
        # this is just an origin translation so without a distance it cannot go ahead
        if isinstance(planet_coo.data, UnitSphericalRepresentation):
            raise u.UnitsError(_NEED_ORIGIN_HINT.format(planet_coo.__class__.__name__))

        if planet_coo.data.differentials:
            bary_sun_pos, bary_sun_vel = get_body_barycentric_posvel(
                planet_coo.body.name, planet_coo.obstime
            )
            bary_sun_pos = bary_sun_pos.with_differentials(
                bary_sun_vel.represent_as(CartesianDifferential)
            )

        else:
            bary_sun_pos = get_body_barycentric(
                planet_coo.body.name, planet_coo.obstime
            )
            bary_sun_vel = None

        return None, bary_sun_pos

    @staticmethod
    def from_icrs(icrs_coo, planet_frame):
        # this is just an origin translation so without a distance it cannot go ahead
        if isinstance(icrs_coo.data, UnitSphericalRepresentation):
            raise u.UnitsError(_NEED_ORIGIN_HINT.format(icrs_coo.__class__.__name__))

        if icrs_coo.data.differentials:
            bary_sun_pos, bary_sun_vel = get_body_barycentric_posvel(
                planet_frame.body.name, planet_frame.obstime
            )
            # Beware! Negation operation is not supported for differentials
            bary_sun_pos = (-bary_sun_pos).with_differentials(
                -bary_sun_vel.represent_as(CartesianDifferential)
            )

        else:
            bary_sun_pos = -get_body_barycentric(
                planet_frame.body.name, planet_frame.obstime
            )
            bary_sun_vel = None

        return None, bary_sun_pos

    @staticmethod
    def self_transform(from_coo, to_frame):
        if np.all(from_coo.obstime == to_frame.obstime):
            return to_frame.realize_frame(from_coo.data)
        else:
            # like CIRS, we do this self-transform via ICRS
            return from_coo.transform_to(ICRS).transform_to(to_frame)


# Redefine HCRS, see https://github.com/astropy/astropy/issues/6835
class HCRS(_PlanetaryICRS):
    body = Sun


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
    )  # negation is required

    # This line works around some input/output quirks of Astropy's rotation_matrix()
    matrix = np.array(rotation_matrix(rotation_angle, rotation_axis.xyz.value.tolist()))
    return matrix


def _obliquity_rotation_value(equinox):
    """
    Function to calculate obliquity of the earth.
    This uses obl06 of erfa.
    """
    jd1, jd2 = get_jd12(equinox, "tt")
    obl = erfa.obl06(jd1, jd2) * u.radian
    return obl.to(u.deg)


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


_FRAME_MAPPING = {
    Sun: {Planes.EARTH_EQUATOR: HCRS, Planes.EARTH_ECLIPTIC: HeliocentricEclipticJ2000},
    Mercury: {Planes.EARTH_EQUATOR: MercuryICRS},
    Venus: {Planes.EARTH_EQUATOR: VenusICRS},
    Earth: {Planes.EARTH_EQUATOR: GCRS},
    Mars: {Planes.EARTH_EQUATOR: MarsICRS},
    Jupiter: {Planes.EARTH_EQUATOR: JupiterICRS},
    Saturn: {Planes.EARTH_EQUATOR: SaturnICRS},
    Uranus: {Planes.EARTH_EQUATOR: UranusICRS},
    Neptune: {Planes.EARTH_EQUATOR: NeptuneICRS},
    Pluto: {Planes.EARTH_EQUATOR: PlutoICRS},
}  # type: Dict[_Body, Dict[Planes, FrameMeta]]


def get_frame(attractor, plane, obstime=J2000):
    """Returns an appropriate reference frame from an attractor and a plane.

    Available planes are Earth equator (parallel to GCRS) and Earth ecliptic.
    The fundamental direction of both is the equinox of epoch (J2000).
    An obstime is needed to properly locate the attractor.

    Parameters
    ----------
    attractor : ~poliastro.bodies.Body
        Body that serves as the center of the frame.
    plane : ~poliastro.frames.Planes
        Fundamental plane of the frame.
    obstime : ~astropy.time.Time
        Time of the frame.

    """
    try:
        frames = _FRAME_MAPPING[attractor]
    except KeyError:
        raise NotImplementedError(
            "Frames for orbits around custom bodies are not yet supported"
        )

    try:
        frame_class = frames[plane]
    except KeyError:
        raise NotImplementedError(
            "A frame with plane {} around body {} is not yet implemented".format(
                plane, attractor
            )
        )

    return frame_class(obstime=obstime)
