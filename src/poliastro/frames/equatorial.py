import numpy as np
from astropy import units as u
from astropy.coordinates import (
    GCRS as _GCRS,
    HCRS as _HCRS,
    ICRS as _ICRS,
    AffineTransform,
    BaseRADecFrame,
    CartesianDifferential,
    FunctionTransformWithFiniteDifference,
    TimeAttribute,
    UnitSphericalRepresentation,
    frame_transform_graph,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME

from poliastro.bodies import (
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Saturn,
    Uranus,
    Venus,
)

__all__ = [
    "ICRS",
    "HCRS",
    "MercuryICRS",
    "VenusICRS",
    "GCRS",
    "MarsICRS",
    "JupiterICRS",
    "SaturnICRS",
    "UranusICRS",
    "NeptuneICRS",
]

# HACK: sphinx-autoapi variable definition
ICRS = _ICRS
HCRS = _HCRS
GCRS = _GCRS


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
        # This is just an origin translation so without a distance it cannot go ahead
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
        # This is just an origin translation so without a distance it cannot go ahead
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
            # Like CIRS, we do this self-transform via ICRS
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


class MoonICRS(_PlanetaryICRS):
    body = Moon


_NEED_ORIGIN_HINT = (
    "The input {0} coordinates do not have length units. This "
    "probably means you created coordinates with lat/lon but "
    "no distance.  PlanetaryICRS<->ICRS transforms cannot "
    "function in this case because there is an origin shift."
)
