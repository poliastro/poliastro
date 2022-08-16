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

Phobos
Deimos


and a way to define new bodies (:py:class:`~Body` class).

Data references can be found in :py:mod:`~poliastro.constants`
"""
import math
from collections import namedtuple

from astropy import units as u
from astropy.constants import G
from astropy.units import Quantity

from poliastro import constants
from poliastro.frames import Planes


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
            "mean_a",
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
        mean_a = 0.0 * u.km,
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
            _q(mean_a),
        )

    @property
    def angular_velocity(self):
        return (2 * math.pi * u.rad) / self.rotational_period.to(u.s)

    def __str__(self):
        return f"{self.name} ({self.symbol})"

    def __reduce__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    @classmethod
    @u.quantity_input(k=u.km**3 / u.s**2, R=u.km)
    def from_parameters(cls, parent, k, name, symbol, R, **kwargs):
        return cls(parent, k, name, symbol, R, **kwargs)

    @classmethod
    def from_relative(
        cls, reference, parent, k, name, symbol=None, R=0, **kwargs
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
        plane : ~poliastro.frames.Planes
            Reference plane of the coordinates.

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
            from poliastro.plotting.interactive import OrbitPlotter3D

            return OrbitPlotter3D(plane=plane).plot_body_orbit(
                self, epoch, label=label
            )
        else:
            from poliastro.plotting.interactive import OrbitPlotter2D

            return OrbitPlotter2D(plane=plane).plot_body_orbit(
                self, epoch, label=label
            )


Sun = Body(
    parent=None,
    k=constants.GM_sun,
    name="Sun",
    symbol="\u2609",
    R=constants.R_sun,
    rotational_period=constants.rotational_period_sun,
    J2=constants.J2_sun,
    mass=constants.M_sun,
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
    mean_a=constants.mean_a_mercury,
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
    J2=constants.J2_venus,
    J3=constants.J3_venus,
    mean_a=constants.mean_a_venus,
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
    mass=constants.M_earth,
    J2=constants.J2_earth,
    J3=constants.J3_earth,
    mean_a=constants.mean_a_earth,
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
    J2=constants.J2_mars,
    J3=constants.J3_mars,
    mean_a=constants.mean_a_mars,
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
    mass=constants.M_jupiter,
    mean_a=constants.mean_a_jupiter,
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
    mean_a=constants.mean_a_saturn,
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
    mean_a=constants.mean_a_uranus,
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
    mean_a=constants.mean_a_neptune,
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
) # No mean_a_pluto as Pluto is officially not a planet around Sun


Moon = Body(
    parent=Earth,
    k=constants.GM_moon,
    name="Moon",
    symbol="\u263E",
    R=constants.R_moon,
    R_mean=constants.R_mean_moon,
    R_polar=constants.R_polar_moon,
    rotational_period=constants.rotational_period_moon,
    mean_a=constants.mean_a_moon
)


Phobos = Body(
    parent=Mars,
    k=constants.GM_phobos,
    name="Phobos",
    mean_a=constants.mean_a_phobos,
)

Deimoms = Body(
    parent=Mars,
    k=constants.GM_deimos,
    name="Deimos",
    mean_a=constants.mean_a_deimos,
)

Europa = Body(
    parent=Jupiter,
    k=constants.GM_europa,
    name="Europa",
    mean_a=constants.mean_a_europa,
)

Ganymede = Body(
    parent=Jupiter,
    k=constants.GM_ganymede,
    name="Ganymede",
    mean_a=constants.mean_a_ganymede,
)

Enceladus = Body(
    parent=Saturn,
    k=constants.GM_enceladus,
    name="Enceladus",
    mean_a=constants.mean_a_enceladus,
)

Titan = Body(
    parent=Saturn,
    k=constants.GM_titan,
    name="Titan",
    mean_a=constants.mean_a_titan,
)

Titania = Body(
    parent=Uranus,
    k=constants.GM_titania,
    name="Titania",
    mean_a=constants.mean_a_titania,
)

Triton = Body(
    parent=Neptune,
    k=constants.GM_triton,
    name="Triton",
    mean_a=constants.mean_a_triton,
)   

Charon = Body(
    parent=Pluto,
    k=constants.GM_charon,
    name="charon",
    mean_a=constants.mean_a_charon,
)






