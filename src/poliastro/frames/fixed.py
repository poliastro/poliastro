from astropy import units as u
from astropy.coordinates import (
    HCRS,
    ITRS as _ITRS,
    BaseRADecFrame,
    FunctionTransform,
    TimeAttribute,
    frame_transform_graph,
)
from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME
from astropy.coordinates.matrix_utilities import rotation_matrix

from poliastro.bodies import (
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Saturn,
    Sun,
    Uranus,
    Venus,
)
from poliastro.constants import J2000
from poliastro.core.fixed import (
    jupiter_rot_elements_at_epoch as jupiter_rot_elements_at_epoch_fast,
    mars_rot_elements_at_epoch as mars_rot_elements_at_epoch_fast,
    mercury_rot_elements_at_epoch as mercury_rot_elements_at_epoch_fast,
    moon_rot_elements_at_epoch as moon_rot_elements_at_epoch_fast,
    neptune_rot_elements_at_epoch as neptune_rot_elements_at_epoch_fast,
    saturn_rot_elements_at_epoch as saturn_rot_elements_at_epoch_fast,
    sun_rot_elements_at_epoch as sun_rot_elements_at_epoch_fast,
    uranus_rot_elements_at_epoch as uranus_rot_elements_at_epoch_fast,
    venus_rot_elements_at_epoch as venus_rot_elements_at_epoch_fast,
)
from poliastro.frames.equatorial import (
    JupiterICRS,
    MarsICRS,
    MercuryICRS,
    MoonICRS,
    NeptuneICRS,
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
]

# HACK: sphinx-autoapi variable definition
ITRS = _ITRS


class _PlanetaryFixed(BaseRADecFrame):
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

    def __new__(cls, *args, **kwargs):
        frame_transform_graph.transform(
            FunctionTransform, cls, cls.equatorial
        )(cls.to_equatorial)
        frame_transform_graph.transform(
            FunctionTransform, cls.equatorial, cls
        )(cls.from_equatorial)

        return super().__new__(cls)

    @staticmethod
    def to_equatorial(fixed_coo, equatorial_frame):
        # TODO replace w/ something smart (Sun/Earth special cased)
        if fixed_coo.body == Sun and type(equatorial_frame) != HCRS:
            raise ValueError(
                f"Equatorial coordinates must be of type `HCRS`, got `{type(equatorial_frame)}` instead."
            )
        elif fixed_coo.body != Sun and fixed_coo.body != equatorial_frame.body:
            raise ValueError(
                "Fixed and equatorial coordinates must have the same body if the fixed frame body is not Sun"
            )

        r = fixed_coo.cartesian

        ra, dec, W = fixed_coo.rot_elements_at_epoch(equatorial_frame.obstime)

        r = r.transform(rotation_matrix(-W, "z"))

        r_trans1 = r.transform(rotation_matrix(-(90 * u.deg - dec), "x"))
        data = r_trans1.transform(rotation_matrix(-(90 * u.deg + ra), "z"))

        return equatorial_frame.realize_frame(data)

    @staticmethod
    def from_equatorial(equatorial_coo, fixed_frame):
        # TODO replace w/ something smart (Sun/Earth special cased)
        if fixed_frame.body == Sun and type(equatorial_coo) != HCRS:
            raise ValueError(
                f"Equatorial coordinates must be of type `HCRS`, got `{type(equatorial_coo)}` instead."
            )
        elif (
            fixed_frame.body != Sun and equatorial_coo.body != fixed_frame.body
        ):
            raise ValueError(
                "Fixed and equatorial coordinates must have the same body if the fixed frame body is not Sun"
            )

        r = equatorial_coo.cartesian

        ra, dec, W = fixed_frame.rot_elements_at_epoch(fixed_frame.obstime)

        r_trans2 = r.transform(rotation_matrix(90 * u.deg + ra, "z"))
        r_f = r_trans2.transform(rotation_matrix(90 * u.deg - dec, "x"))
        r_f = r_f.transform(rotation_matrix(W, "z"))

        return fixed_frame.realize_frame(r_f)

    @classmethod
    def rot_elements_at_epoch(cls, epoch=J2000):
        """Provides rotational elements at epoch.

        Provides north pole of body and angle to prime meridian.

        Parameters
        ----------
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        Returns
        -------
        ra, dec, W: tuple (~astropy.units.Quantity)
            Right ascension and declination of north pole, and angle of the prime meridian.

        """
        T = (epoch.tdb - J2000).to_value(u.d) / 36525
        d = (epoch.tdb - J2000).to_value(u.d)
        return cls._rot_elements_at_epoch(T, d)

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        raise NotImplementedError


class SunFixed(_PlanetaryFixed):
    body = Sun
    equatorial = HCRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = sun_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class MercuryFixed(_PlanetaryFixed):
    body = Mercury
    equatorial = MercuryICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = mercury_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class VenusFixed(_PlanetaryFixed):
    body = Venus
    equatorial = VenusICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = venus_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class MarsFixed(_PlanetaryFixed):
    body = Mars
    equatorial = MarsICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = mars_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class JupiterFixed(_PlanetaryFixed):
    body = Jupiter
    equatorial = JupiterICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = jupiter_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class SaturnFixed(_PlanetaryFixed):
    body = Saturn
    equatorial = SaturnICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = saturn_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class UranusFixed(_PlanetaryFixed):
    body = Uranus
    equatorial = UranusICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = uranus_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class NeptuneFixed(_PlanetaryFixed):
    body = Neptune
    equatorial = NeptuneICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = neptune_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg


class MoonFixed(_PlanetaryFixed):
    body = Moon
    equatorial = MoonICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra, dec, W = moon_rot_elements_at_epoch_fast(T, d)
        return ra * u.deg, dec * u.deg, W * u.deg
