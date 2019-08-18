from astropy.coordinates import (
    HCRS,
    ITRS,
    BaseRADecFrame,
    CartesianDifferential,
    CartesianRepresentation,
    FunctionTransform,
    TimeAttribute,
    frame_transform_graph,
)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME

from poliastro.bodies import (
    Jupiter,
    Mars,
    Mercury,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
)
from poliastro.coordinates import body_centered_to_icrs, icrs_to_body_centered

from .equatorial import (
    JupiterICRS,
    MarsICRS,
    MercuryICRS,
    NeptuneICRS,
    PlutoICRS,
    SaturnICRS,
    UranusICRS,
    VenusICRS,
)

__all__ = [
    "SunFixed",
    "MercuryFixed",
    "VenusFixed",
    "ITRS",
    "MarsFixed",
    "JupiterFixed",
    "SaturnFixed",
    "UranusFixed",
    "NeptuneFixed",
    "PlutoFixed",
]


class _PlanetaryFixed(BaseRADecFrame):
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(FunctionTransform, cls, cls.equatorial)(
            cls.to_equatorial
        )
        frame_transform_graph.transform(FunctionTransform, cls.equatorial, cls)(
            cls.from_equatorial
        )

        return super().__new__(cls)

    @staticmethod
    def to_equatorial(fixed_coo, equatorial_frame):
        assert fixed_coo.body == equatorial_frame.body

        r, v = body_centered_to_icrs(
            fixed_coo.data.xyz,
            fixed_coo.data.differentials["s"].d_xyz,
            fixed_coo.body,
            fixed_coo.obstime,
            rotate_meridian=True,
        )
        data = CartesianRepresentation(r, differentials=CartesianDifferential(v))
        return equatorial_frame.realize_frame(data)

    @staticmethod
    def from_equatorial(equatorial_coo, fixed_frame):
        assert equatorial_coo.body == fixed_frame.body

        r, v = icrs_to_body_centered(
            equatorial_coo.data.xyz,
            equatorial_coo.data.differentials["s"].d_xyz,
            equatorial_coo.body,
            equatorial_coo.obstime,
            rotate_meridian=True,
        )
        data = CartesianRepresentation(r, differentials=CartesianDifferential(v))
        return fixed_frame.realize_frame(data)

class SunFixed(_PlanetaryFixed):
    body = Sun
    equatorial = HCRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 286.13 * u.deg
        dec = 63.87 * u.deg
        W = (84.176 + 14.1844000 * d) * u.deg

        return ra, dec, W


class MercuryFixed(_PlanetaryFixed):
    body = Mercury
    equatorial = MercuryICRS


class VenusFixed(_PlanetaryFixed):
    body = Venus
    equatorial = VenusICRS


class MarsFixed(_PlanetaryFixed):
    body = Mars
    equatorial = MarsICRS


class JupiterFixed(_PlanetaryFixed):
    body = Jupiter
    equatorial = JupiterICRS


class SaturnFixed(_PlanetaryFixed):
    body = Saturn
    equatorial = SaturnICRS


class UranusFixed(_PlanetaryFixed):
    body = Uranus
    equatorial = UranusICRS


class NeptuneFixed(_PlanetaryFixed):
    body = Neptune
    equatorial = NeptuneICRS


class PlutoFixed(_PlanetaryFixed):
    body = Pluto
    equatorial = PlutoICRS
