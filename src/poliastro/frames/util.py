from typing import Dict

from astropy.coordinates.baseframe import FrameMeta

from poliastro.bodies import (
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Neptune,
    Pluto,
    Saturn,
    SolarSystemBody,
    Sun,
    Uranus,
    Venus,
)
from poliastro.constants import J2000

from .ecliptic import GeocentricMeanEcliptic, HeliocentricEclipticJ2000
from .enums import Planes
from .equatorial import (
    GCRS,
    HCRS,
    JupiterICRS,
    MarsICRS,
    MercuryICRS,
    NeptuneICRS,
    PlutoICRS,
    SaturnICRS,
    UranusICRS,
    VenusICRS,
)

_FRAME_MAPPING = {
    Sun: {Planes.EARTH_EQUATOR: HCRS, Planes.EARTH_ECLIPTIC: HeliocentricEclipticJ2000},
    Mercury: {Planes.EARTH_EQUATOR: MercuryICRS},
    Venus: {Planes.EARTH_EQUATOR: VenusICRS},
    Earth: {Planes.EARTH_EQUATOR: GCRS, Planes.EARTH_ECLIPTIC: GeocentricMeanEcliptic},
    Mars: {Planes.EARTH_EQUATOR: MarsICRS},
    Jupiter: {Planes.EARTH_EQUATOR: JupiterICRS},
    Saturn: {Planes.EARTH_EQUATOR: SaturnICRS},
    Uranus: {Planes.EARTH_EQUATOR: UranusICRS},
    Neptune: {Planes.EARTH_EQUATOR: NeptuneICRS},
    Pluto: {Planes.EARTH_EQUATOR: PlutoICRS},
}  # type: Dict[SolarSystemBody, Dict[Planes, FrameMeta]]


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
