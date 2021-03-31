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
from collections import namedtuple

from astropy import units as u
from astropy.constants import G
from astropy.units import Quantity

from . import constants
from .frames import Planes


# HACK: Constants cannot be hashed
# (see https://github.com/astropy/astropy/issues/10043)
# so we will convert them all to normal Quantities
def _q(c):
    return Quantity(c)


# https://stackoverflow.com/a/16721002/554319
class Body(
    namedtuple(
        "_Body",
        [
            "parent",
            "k",
            "name",
            "symbol",
            "R",
            "R_polar",
            "R_mean",
            "rotational_period",
            "J2",
            "J3",
            "mass",
        ],
    )
):
    __slots__ = ()

    def __new__(
        cls,
        parent,
        k,
        name,
        symbol=None,
        R=0 * u.km,
        R_polar=0 * u.km,
        R_mean=0 * u.km,
        rotational_period=0.0 * u.day,
        J2=0.0 * u.one,
        J3=0.0 * u.one,
        mass=None,
    ):
        if mass is None:
            mass = k / G

        return super().__new__(
            cls,
            parent,
            _q(k),
            name,
            symbol,
            _q(R),
            _q(R_polar),
            _q(R_mean),
            _q(rotational_period),
            _q(J2),
            _q(J3),
            _q(mass),
        )

    @property
    def angular_velocity(self):
        return (2 * math.pi * u.rad) / self.rotational_period.to(u.s)

    def __str__(self):
        return f"{self.name} ({self.symbol})"

    def __repr__(self):
        return self.__str__()

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


class SolarSystemPlanet(Body):
    def plot(
        self,
        epoch=None,
        label=None,
        use_3d=False,
        interactive=False,
        plane=Planes.EARTH_ECLIPTIC,
    ):
        """Plots the body orbit.

        Parameters
        ----------
        epoch : astropy.time.Time, optional
            Epoch of current position.
        label : str, optional
            Label for the orbit, defaults to empty.
        use_3d : bool, optional
            Produce a 3D plot, default to False.
        interactive : bool, optional
            Produce an interactive (rather than static) image of the orbit, default to False.
            This option requires Plotly properly installed and configured for your environment.

        """
        if not interactive and use_3d:
            raise ValueError(
                "The static plotter does not support 3D, use `interactive=True`"
            )
        elif not interactive:
            from poliastro.plotting.static import StaticOrbitPlotter

            return StaticOrbitPlotter(plane=plane).plot_body_orbit(
                self, epoch, label=label
            )
        elif use_3d:
            from poliastro.plotting.core import OrbitPlotter3D

            return OrbitPlotter3D(plane=plane).plot_body_orbit(self, epoch, label=label)
        else:
            from poliastro.plotting.core import OrbitPlotter2D

            return OrbitPlotter2D(plane=plane).plot_body_orbit(self, epoch, label=label)


Sun = SolarSystemPlanet(
    parent=None,
    k=constants.GM_sun,
    name="Sun",
    symbol="\u2609",
    R=constants.R_sun,
    rotational_period=constants.rotational_period_sun,
    J2=_q(constants.J2_sun),
    mass=_q(constants.M_sun),
)


Mercury = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_mercury,
    name="Mercury",
    symbol="\u263F",
    R=constants.R_mercury,
    R_mean=constants.R_mean_mercury,
    R_polar=constants.R_polar_mercury,
    rotational_period=constants.rotational_period_mercury,
)

Venus = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_venus,
    name="Venus",
    symbol="\u2640",
    R=constants.R_venus,
    R_mean=constants.R_mean_venus,
    R_polar=constants.R_polar_venus,
    rotational_period=constants.rotational_period_venus,
    J2=_q(constants.J2_venus),
    J3=_q(constants.J3_venus),
)
Earth = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_earth,
    name="Earth",
    symbol="\u2641",
    R=constants.R_earth,
    R_mean=constants.R_mean_earth,
    R_polar=constants.R_polar_earth,
    rotational_period=constants.rotational_period_earth,
    mass=_q(constants.M_earth),
    J2=_q(constants.J2_earth),
    J3=_q(constants.J3_earth),
)
Mars = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_mars,
    name="Mars",
    symbol="\u2642",
    R=constants.R_mars,
    R_mean=constants.R_mean_mars,
    R_polar=constants.R_polar_mars,
    rotational_period=constants.rotational_period_mars,
    J2=_q(constants.J2_mars),
    J3=_q(constants.J3_mars),
)
Jupiter = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_jupiter,
    name="Jupiter",
    symbol="\u2643",
    R=constants.R_jupiter,
    R_mean=constants.R_mean_jupiter,
    R_polar=constants.R_polar_jupiter,
    rotational_period=constants.rotational_period_jupiter,
    mass=_q(constants.M_jupiter),
)
Saturn = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_saturn,
    name="Saturn",
    symbol="\u2644",
    R=constants.R_saturn,
    R_mean=constants.R_mean_saturn,
    R_polar=constants.R_polar_saturn,
    rotational_period=constants.rotational_period_saturn,
)
Uranus = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_uranus,
    name="Uranus",
    symbol="\u26E2",
    R=constants.R_uranus,
    R_mean=constants.R_mean_uranus,
    R_polar=constants.R_polar_uranus,
    rotational_period=constants.rotational_period_uranus,
)
Neptune = SolarSystemPlanet(
    parent=Sun,
    k=constants.GM_neptune,
    name="Neptune",
    symbol="\u2646",
    R=constants.R_neptune,
    R_mean=constants.R_mean_neptune,
    R_polar=constants.R_polar_neptune,
    rotational_period=constants.rotational_period_neptune,
)

Pluto = Body(
    parent=Sun,
    k=constants.GM_pluto,
    name="Pluto",
    symbol="\u2647",
    R=constants.R_pluto,
    R_mean=constants.R_mean_pluto,
    R_polar=constants.R_polar_pluto,
    rotational_period=constants.rotational_period_pluto,
)

Moon = Body(
    parent=Earth,
    k=constants.GM_moon,
    name="Moon",
    symbol="\u263E",
    R=constants.R_moon,
    R_mean=constants.R_mean_moon,
    R_polar=constants.R_polar_moon,
    rotational_period=constants.rotational_period_moon,
)
