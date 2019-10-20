"""Bodies of the Solar System.

Contains some predefined bodies of the Solar System:

* Sun (☉)
* Earth (♁)
* Moon (☾)
* Mercury (☿)
* Venus (♀)
* Mars (♂)
* Jupiter (♃)
* Saturn (♄)
* Uranus (⛢)
* Neptune (♆)
* Pluto (♇)

and a way to define new bodies (:py:class:`~Body` class).

Data references can be found in :py:mod:`~poliastro.constants`
"""
import math

from astropy import units as u
from astropy.constants import G

from poliastro import constants


class _Body(object):
    name: str
    symbol: str

    def __str__(self):
        return f"{self.name} ({self.symbol})"

    def __repr__(self):
        return self.__str__()

    def rot_elements_at_epoch(self, epoch):
        """Provides rotational elements at epoch.

        Provides north pole of body and angle to prime meridian

        Parameters
        ----------
        epoch : ~astropy.time.Time, optional
            Epoch, default to J2000.

        Returns
        -------
        ra, dec, W: tuple (~astropy.units.Quantity)
            Right ascension and declination of north pole, and angle of the prime meridian.

        """
        T = (epoch.tdb - constants.J2000).to("day").value / 36525
        d = (epoch.tdb - constants.J2000).to("day").value
        return self._rot_elements_at_epoch(T, d)

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        raise NotImplementedError("Function only defined for some Solar System bodies")


class Body(_Body):
    """Class to represent a generic body.

    """

    def __init__(
        self,
        parent,
        k,
        name,
        symbol=None,
        R=0 * u.km,
        R_polar=0 * u.km,
        R_mean=0 * u.km,
        rotational_period=None,
        **kwargs,
    ):
        """Constructor.

        Parameters
        ----------
        parent : Body
            Central body.
        k : ~astropy.units.Quantity
            Standard gravitational parameter.
        name : str
            Name of the body.
        symbol : str, optional
            Symbol for the body.
        R : ~astropy.units.Quantity, optional
            Equatorial radius of the body.
        R_polar: ~astropy.units.Quantity, optional
            Polar radius of the body
        R_mean: ~astropy.units.Quantity, optional
            Mean radius of the body
        rotational_period: ~astropy.units.Quantity, optional
            Rotational period
        """
        self.parent = parent
        self.k = k
        self.name = name
        self.symbol = symbol
        self.R = R
        self.R_polar = R_polar
        self.R_mean = R_mean
        self.rotational_period = rotational_period
        self.kwargs = kwargs

    @property
    def angular_velocity(self):
        return (2 * math.pi * u.rad) / (self.rotational_period).to(u.s)

    @classmethod
    @u.quantity_input(k=u.km ** 3 / u.s ** 2, R=u.km)
    def from_parameters(cls, parent, k, name, symbol, R, **kwargs):
        return cls(parent, k, name, symbol, R, **kwargs)

    @classmethod
    def from_relative(
        cls, reference, parent=None, k=None, name=None, symbol=None, R=None, **kwargs
    ):
        k = k * reference.k
        R = R * reference.R
        return cls(parent, k, name, symbol, R, **kwargs)


class _Sun(_Body):
    parent = None
    k = constants.GM_sun
    name = "Sun"
    symbol = "\u2609"
    R = constants.R_sun
    mass = k / G
    J2 = constants.J2_sun
    rotational_period = constants.rotational_period_sun
    Wdivc = constants.Wdivc_sun

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 286.13 * u.deg
        dec = 63.87 * u.deg
        W = (84.176 + 14.1844000 * d) * u.deg

        return ra, dec, W


Sun = _Sun()


class _Mercury(_Body):
    parent = Sun
    k = constants.GM_mercury
    name = "Mercury"
    symbol = "\u263F"
    R = constants.R_mercury
    R_mean = constants.R_mean_mercury
    R_polar = constants.R_polar_mercury
    rotational_period = constants.rotational_period_mercury
    mass = k / G

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
            0.01067257 * math.sin(M1.to("rad").value)
            - 0.00112309 * math.sin(M2.to("rad").value)
            - 0.00011040 * math.sin(M3.to("rad").value)
            - 0.00002539 * math.sin(M4.to("rad").value)
            - 0.00000571 * math.sin(M5.to("rad").value)
        ) * u.deg

        return ra, dec, W


class _Venus(_Body):
    parent = Sun
    k = constants.GM_venus
    name = "Venus"
    symbol = "\u2640"
    R = constants.R_venus
    R_mean = constants.R_mean_venus
    R_polar = constants.R_polar_venus
    rotational_period = constants.rotational_period_venus
    mass = k / G
    J2 = constants.J2_venus
    J3 = constants.J3_venus
    ecc = 0.007 * u.one

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 272.76 * u.deg
        dec = 67.16 * u.deg
        W = (160.20 - 1.4813688 * d) * u.deg

        return ra, dec, W


class _Earth(_Body):
    parent = Sun
    k = constants.GM_earth
    name = "Earth"
    symbol = "\u2641"
    R = constants.R_earth
    R_mean = constants.R_mean_earth
    R_polar = constants.R_polar_earth
    rotational_period = constants.rotational_period_earth
    mass = k / G
    J2 = constants.J2_earth
    J3 = constants.J3_earth
    H0 = constants.H0_earth
    rho0 = constants.rho0_earth

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = (0.00 - 0.641 * T) * u.deg
        dec = (90.00 - 0.557 * T) * u.deg
        W = (190.147 + 360.9856235 * d) * u.deg

        return ra, dec, W


class _Mars(_Body):
    parent = Sun
    k = constants.GM_mars
    name = "Mars"
    symbol = "\u2642"
    R = constants.R_mars
    R_mean = constants.R_mean_mars
    R_polar = constants.R_polar_mars
    rotational_period = constants.rotational_period_mars
    mass = k / G
    J2 = constants.J2_mars
    J3 = constants.J3_mars

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
            + 0.000068 * math.sin(M1.to("rad").value)
            + 0.000238 * math.sin(M2.to("rad").value)
            + 0.000052 * math.sin(M3.to("rad").value)
            + 0.000009 * math.sin(M4.to("rad").value)
            + 0.419057 * math.sin(M5.to("rad").value)
        ) * u.deg

        K1 = (122.433576 + 19139.9407476 * T) * u.deg
        K2 = (43.058401 + 38280.8753272 * T) * u.deg
        K3 = (57.663379 + 57420.7517205 * T) * u.deg
        K4 = (79.476401 + 76560.6495004 * T) * u.deg
        K5 = (166.325722 + 0.5042615 * T) * u.deg

        dec = (
            54.432516
            - 0.05827105 * T
            + 0.000051 * math.cos(K1.to("rad").value)
            + 0.000141 * math.cos(K2.to("rad").value)
            + 0.000031 * math.cos(K3.to("rad").value)
            + 0.000005 * math.cos(K4.to("rad").value)
            + 1.591274 * math.cos(K5.to("rad").value)
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
            + 0.000145 * math.sin(J1.to("rad").value)
            + 0.000157 * math.sin(J2.to("rad").value)
            + 0.000040 * math.sin(J3.to("rad").value)
            + 0.000001 * math.sin(J4.to("rad").value)
            + 0.000001 * math.sin(J5.to("rad").value)
            + 0.584542 * math.sin(J6.to("rad").value)
        ) * u.deg

        return ra, dec, W


class _Jupiter(_Body):
    parent = Sun
    k = constants.GM_jupiter
    name = "Jupiter"
    symbol = "\u2643"
    R = constants.R_jupiter
    R_mean = constants.R_mean_jupiter
    R_polar = constants.R_polar_jupiter
    rotational_period = constants.rotational_period_jupiter
    mass = k / G

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
            + 0.000117 * math.sin(Ja.to("rad").value)
            + 0.000938 * math.sin(Jb.to("rad").value)
            + 0.001432 * math.sin(Jc.to("rad").value)
            + 0.000030 * math.sin(Jd.to("rad").value)
            + 0.002150 * math.sin(Je.to("rad").value)
        ) * u.deg
        dec = (
            64.495303
            + 0.002413 * T
            + 0.000050 * math.cos(Ja.to("rad").value)
            + 0.000404 * math.cos(Jb.to("rad").value)
            + 0.000617 * math.cos(Jc.to("rad").value)
            - 0.000013 * math.cos(Jd.to("rad").value)
            + 0.000926 * math.cos(Je.to("rad").value)
        ) * u.deg
        W = (284.95 + 870.536 * d) * u.deg

        return ra, dec, W


class _Saturn(_Body):
    parent = Sun
    k = constants.GM_saturn
    name = "Saturn"
    symbol = "\u2644"
    R = constants.R_saturn
    R_mean = constants.R_mean_saturn
    R_polar = constants.R_polar_saturn
    rotational_period = constants.rotational_period_saturn
    mass = k / G

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = (40.589 - 0.036 * T) * u.deg
        dec = (83.537 - 0.004 * T) * u.deg
        W = (38.90 + 810.7939024 * d) * u.deg

        return ra, dec, W


class _Uranus(_Body):
    parent = Sun
    k = constants.GM_uranus
    name = "Uranus"
    symbol = "\u26E2"
    R = constants.R_uranus
    R_mean = constants.R_mean_uranus
    R_polar = constants.R_polar_uranus
    rotational_period = constants.rotational_period_uranus
    mass = k / G

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 257.311 * u.deg
        dec = -15.175 * u.deg
        W = (203.81 - 501.1600928 * d) * u.deg

        return ra, dec, W


class _Neptune(_Body):
    parent = Sun
    k = constants.GM_neptune
    name = "Neptune"
    symbol = "\u2646"
    R = constants.R_neptune
    R_mean = constants.R_mean_neptune
    R_polar = constants.R_polar_neptune
    rotational_period = constants.rotational_period_neptune
    mass = k / G

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        N = (357.85 + 52.316 * T) * u.deg

        ra = (299.36 + 0.70 * math.sin(N.to("rad").value)) * u.deg
        dec = (43.46 - 0.51 * math.cos(N.to("rad").value)) * u.deg
        W = (249.978 + 541.1397757 * d - 0.48 * math.sin(N.to("rad").value)) * u.deg

        return ra, dec, W


class _Pluto(_Body):
    parent = Sun
    k = constants.GM_pluto
    name = "Pluto"
    symbol = "\u2647"
    R = constants.R_pluto
    R_mean = constants.R_mean_pluto
    R_polar = constants.R_polar_pluto
    rotational_period = constants.rotational_period_pluto
    mass = k / G

    @staticmethod
    def _rot_elements_at_epoch(T, d):
        ra = 132.993 * u.deg
        dec = -6.163 * u.deg
        W = (302.695 + 56.3625225 * d) * u.deg

        return ra, dec, W


Mercury = _Mercury()
Venus = _Venus()
Earth = _Earth()
Mars = _Mars()
Jupiter = _Jupiter()
Saturn = _Saturn()
Uranus = _Uranus()
Neptune = _Neptune()

Pluto = _Pluto()


class _Moon(_Body):
    parent = Earth
    k = constants.GM_moon
    name = "Moon"
    symbol = "\u263E"
    R = constants.R_moon
    R_mean = constants.R_mean_moon
    R_polar = constants.R_polar_moon
    rotational_period = constants.rotational_period_moon
    mass = k / G

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
            - 3.8787 * math.sin(E1.to("rad").value)
            - 0.1204 * math.sin(E2.to("rad").value)
            + 0.0700 * math.sin(E3.to("rad").value)
            - 0.0172 * math.sin(E4.to("rad").value)
            + 0.0072 * math.sin(E6.to("rad").value)
            - 0.0052 * math.sin(E10.to("rad").value)
            + 0.0043 * math.sin(E13.to("rad").value)
        ) * u.deg
        dec = (
            66.5392
            + 0.0130 * T
            + 1.5419 * math.cos(E1.to("rad").value)
            + 0.0239 * math.cos(E2.to("rad").value)
            - 0.0278 * math.cos(E3.to("rad").value)
            + 0.0068 * math.cos(E4.to("rad").value)
            - 0.0029 * math.cos(E6.to("rad").value)
            + 0.0009 * math.cos(E7.to("rad").value)
            + 0.0008 * math.cos(E10.to("rad").value)
            - 0.0009 * math.cos(E13.to("rad").value)
        ) * u.deg
        W = (
            38.321
            + 13.17635815 * d
            - 1.4e-12 * d ** 2
            + 3.5610 * math.sin(E1.to("rad").value)
            + 0.1208 * math.sin(E2.to("rad").value)
            - 0.0642 * math.sin(E3.to("rad").value)
            + 0.0158 * math.sin(E4.to("rad").value)
            + 0.0252 * math.sin(E5.to("rad").value)
            - 0.0066 * math.sin(E6.to("rad").value)
            - 0.0047 * math.sin(E7.to("rad").value)
            - 0.0046 * math.sin(E8.to("rad").value)
            + 0.0028 * math.sin(E9.to("rad").value)
            + 0.0052 * math.sin(E10.to("rad").value)
            + 0.0040 * math.sin(E11.to("rad").value)
            + 0.0019 * math.sin(E12.to("rad").value)
            - 0.0044 * math.sin(E13.to("rad").value)
        ) * u.deg

        return ra, dec, W


Moon = _Moon()
