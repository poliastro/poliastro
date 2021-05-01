import numpy as np
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

from .equatorial import (
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
        frame_transform_graph.transform(FunctionTransform, cls, cls.equatorial)(
            cls.to_equatorial
        )
        frame_transform_graph.transform(FunctionTransform, cls.equatorial, cls)(
            cls.from_equatorial
        )

        return super().__new__(cls)

    @staticmethod
    def to_equatorial(fixed_coo, equatorial_frame):
        # TODO replace w/ something smart (Sun/Earth special cased)
        if fixed_coo.body == Sun:
            assert type(equatorial_frame) == HCRS
        else:
            assert fixed_coo.body == equatorial_frame.body

        r = fixed_coo.cartesian

        ra, dec, W = fixed_coo.rot_elements_at_epoch(equatorial_frame.obstime)

        r = r.transform(rotation_matrix(-W, "z"))

        r_trans1 = r.transform(rotation_matrix(-(90 * u.deg - dec), "x"))
        data = r_trans1.transform(rotation_matrix(-(90 * u.deg + ra), "z"))

        return equatorial_frame.realize_frame(data)

    @staticmethod
    def from_equatorial(equatorial_coo, fixed_frame):
        # TODO replace w/ something smart (Sun/Earth special cased)
        if fixed_frame.body == Sun:
            assert type(equatorial_coo) == HCRS
        else:
            assert equatorial_coo.body == fixed_frame.body

        r = equatorial_coo.cartesian

        ra, dec, W = fixed_frame.rot_elements_at_epoch(fixed_frame.obstime)

        r_trans2 = r.transform(rotation_matrix(90 * u.deg + ra, "z"))
        r_f = r_trans2.transform(rotation_matrix(90 * u.deg - dec, "x"))
        r_f = r_f.transform(rotation_matrix(W, "z"))

        return fixed_frame.realize_frame(r_f)

    @classmethod
    def rot_elements_at_epoch(cls, epoch):
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
        T = (epoch.tdb - J2000).to(u.day).value / 36525
        d = (epoch.tdb - J2000).to(u.day).value
        return cls._rot_elements_at_epoch(T, d)

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        raise NotImplementedError


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

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        M1 = (174.7910857 + 4.092335 * d) * u.deg
        M2 = (349.5821714 + 8.184670 * d) * u.deg
        M3 = (164.3732571 + 12.277005 * d) * u.deg
        M4 = (339.1643429 + 16.369340 * d) * u.deg
        M5 = (153.9554286 + 20.461675 * d) * u.deg
        ra = (281.0103 - 0.0328 * T) * u.deg
        dec = (61.45 - 0.005 * T) * u.deg
        W = (329.5988 + 6.1385108 * d) * u.deg + (
            0.01067257 * np.sin(M1.to("rad").value)
            - 0.00112309 * np.sin(M2.to("rad").value)
            - 0.00011040 * np.sin(M3.to("rad").value)
            - 0.00002539 * np.sin(M4.to("rad").value)
            - 0.00000571 * np.sin(M5.to("rad").value)
        ) * u.deg

        return ra, dec, W


class VenusFixed(_PlanetaryFixed):
    body = Venus
    equatorial = VenusICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 272.76 * u.deg
        dec = 67.16 * u.deg
        W = (160.20 - 1.4813688 * d) * u.deg

        return ra, dec, W


class MarsFixed(_PlanetaryFixed):
    body = Mars
    equatorial = MarsICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        M1 = (198.991226 + 19139.4819985 * T) * u.deg
        M2 = (226.292679 + 38280.8511281 * T) * u.deg
        M3 = (249.663391 + 57420.7251593 * T) * u.deg
        M4 = (266.183510 + 76560.6367950 * T) * u.deg
        M5 = (79.398797 + 0.5042615 * T) * u.deg

        ra = (
            317.269202
            - 0.10927547 * T
            + 0.000068 * np.sin(M1.to("rad").value)
            + 0.000238 * np.sin(M2.to("rad").value)
            + 0.000052 * np.sin(M3.to("rad").value)
            + 0.000009 * np.sin(M4.to("rad").value)
            + 0.419057 * np.sin(M5.to("rad").value)
        ) * u.deg

        K1 = (122.433576 + 19139.9407476 * T) * u.deg
        K2 = (43.058401 + 38280.8753272 * T) * u.deg
        K3 = (57.663379 + 57420.7517205 * T) * u.deg
        K4 = (79.476401 + 76560.6495004 * T) * u.deg
        K5 = (166.325722 + 0.5042615 * T) * u.deg

        dec = (
            54.432516
            - 0.05827105 * T
            + 0.000051 * np.cos(K1.to("rad").value)
            + 0.000141 * np.cos(K2.to("rad").value)
            + 0.000031 * np.cos(K3.to("rad").value)
            + 0.000005 * np.cos(K4.to("rad").value)
            + 1.591274 * np.cos(K5.to("rad").value)
        ) * u.deg

        J1 = (129.071773 + 19140.0328244 * T) * u.deg
        J2 = (36.352167 + 38281.0473591 * T) * u.deg
        J3 = (56.668646 + 57420.9295360 * T) * u.deg
        J4 = (67.364003 + 76560.2552215 * T) * u.deg
        J5 = (104.792680 + 95700.4387578 * T) * u.deg
        J6 = (95.391654 + 0.5042615 * T) * u.deg

        W = (
            176.049863
            + 350.891982443297 * d
            + 0.000145 * np.sin(J1.to("rad").value)
            + 0.000157 * np.sin(J2.to("rad").value)
            + 0.000040 * np.sin(J3.to("rad").value)
            + 0.000001 * np.sin(J4.to("rad").value)
            + 0.000001 * np.sin(J5.to("rad").value)
            + 0.584542 * np.sin(J6.to("rad").value)
        ) * u.deg

        return ra, dec, W


class JupiterFixed(_PlanetaryFixed):
    body = Jupiter
    equatorial = JupiterICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        Ja = (99.360714 + 4850.4046 * T) * u.deg
        Jb = (175.895369 + 1191.9605 * T) * u.deg
        Jc = (300.323162 + 262.5475 * T) * u.deg
        Jd = (114.012305 + 6070.2476 * T) * u.deg
        Je = (49.511251 + 64.3000 * T) * u.deg

        ra = (
            268.056595
            - 0.006499 * T
            + 0.000117 * np.sin(Ja.to("rad").value)
            + 0.000938 * np.sin(Jb.to("rad").value)
            + 0.001432 * np.sin(Jc.to("rad").value)
            + 0.000030 * np.sin(Jd.to("rad").value)
            + 0.002150 * np.sin(Je.to("rad").value)
        ) * u.deg
        dec = (
            64.495303
            + 0.002413 * T
            + 0.000050 * np.cos(Ja.to("rad").value)
            + 0.000404 * np.cos(Jb.to("rad").value)
            + 0.000617 * np.cos(Jc.to("rad").value)
            - 0.000013 * np.cos(Jd.to("rad").value)
            + 0.000926 * np.cos(Je.to("rad").value)
        ) * u.deg
        W = (284.95 + 870.536 * d) * u.deg

        return ra, dec, W


class SaturnFixed(_PlanetaryFixed):
    body = Saturn
    equatorial = SaturnICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = (40.589 - 0.036 * T) * u.deg
        dec = (83.537 - 0.004 * T) * u.deg
        W = (38.90 + 810.7939024 * d) * u.deg

        return ra, dec, W


class UranusFixed(_PlanetaryFixed):
    body = Uranus
    equatorial = UranusICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 257.311 * u.deg
        dec = -15.175 * u.deg
        W = (203.81 - 501.1600928 * d) * u.deg

        return ra, dec, W


class NeptuneFixed(_PlanetaryFixed):
    body = Neptune
    equatorial = NeptuneICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        N = (357.85 + 52.316 * T) * u.deg

        ra = (299.36 + 0.70 * np.sin(N.to("rad").value)) * u.deg
        dec = (43.46 - 0.51 * np.cos(N.to("rad").value)) * u.deg
        W = (249.978 + 541.1397757 * d - 0.48 * np.sin(N.to("rad").value)) * u.deg

        return ra, dec, W


class MoonFixed(_PlanetaryFixed):
    body = Moon
    equatorial = MoonICRS

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        E1 = (125.045 - 0.0529921 * d) * u.deg
        E2 = (250.089 - 0.1059842 * d) * u.deg
        E3 = (260.008 + 13.0120009 * d) * u.deg
        E4 = (176.625 + 13.3407154 * d) * u.deg
        E5 = (357.529 + 0.9856003 * d) * u.deg
        E6 = (311.589 + 26.4057084 * d) * u.deg
        E7 = (134.963 + 13.0649930 * d) * u.deg
        E8 = (276.617 + 0.3287146 * d) * u.deg
        E9 = (34.226 + 1.7484877 * d) * u.deg
        E10 = (15.134 - 0.1589763 * d) * u.deg
        E11 = (119.743 + 0.0036096 * d) * u.deg
        E12 = (239.961 + 0.1643573 * d) * u.deg
        E13 = (25.053 + 12.9590088 * d) * u.deg

        ra = (
            269.9949
            + 0.0031 * T
            - 3.8787 * np.sin(E1.to("rad").value)
            - 0.1204 * np.sin(E2.to("rad").value)
            + 0.0700 * np.sin(E3.to("rad").value)
            - 0.0172 * np.sin(E4.to("rad").value)
            + 0.0072 * np.sin(E6.to("rad").value)
            - 0.0052 * np.sin(E10.to("rad").value)
            + 0.0043 * np.sin(E13.to("rad").value)
        ) * u.deg
        dec = (
            66.5392
            + 0.0130 * T
            + 1.5419 * np.cos(E1.to("rad").value)
            + 0.0239 * np.cos(E2.to("rad").value)
            - 0.0278 * np.cos(E3.to("rad").value)
            + 0.0068 * np.cos(E4.to("rad").value)
            - 0.0029 * np.cos(E6.to("rad").value)
            + 0.0009 * np.cos(E7.to("rad").value)
            + 0.0008 * np.cos(E10.to("rad").value)
            - 0.0009 * np.cos(E13.to("rad").value)
        ) * u.deg
        W = (
            38.3213
            + 13.17635815 * d
            - 1.4e-12 * d ** 2
            + 3.5610 * np.sin(E1.to("rad").value)
            + 0.1208 * np.sin(E2.to("rad").value)
            - 0.0642 * np.sin(E3.to("rad").value)
            + 0.0158 * np.sin(E4.to("rad").value)
            + 0.0252 * np.sin(E5.to("rad").value)
            - 0.0066 * np.sin(E6.to("rad").value)
            - 0.0047 * np.sin(E7.to("rad").value)
            - 0.0046 * np.sin(E8.to("rad").value)
            + 0.0028 * np.sin(E9.to("rad").value)
            + 0.0052 * np.sin(E10.to("rad").value)
            + 0.0040 * np.sin(E11.to("rad").value)
            + 0.0019 * np.sin(E12.to("rad").value)
            - 0.0044 * np.sin(E13.to("rad").value)
        ) * u.deg

        return ra, dec, W
