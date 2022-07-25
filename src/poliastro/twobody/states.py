from abc import ABC, abstractmethod
from functools import cached_property

from astropy import units as u

from poliastro.core.elements import coe2mee, coe2rv, mee2coe, mee2rv, rv2coe
from poliastro.twobody.elements import mean_motion, period, t_p


class BaseState(ABC):
    """Base State class, meant to be subclassed."""

    def __init__(self, attractor, elements, plane, _units=True):
        """Constructor.

        Parameters
        ----------
        attractor : Body
            Main attractor.
        elements : tuple
            Six-tuple of orbital elements for this state.
        plane : ~poliastro.frames.enums.Planes
            Reference plane for the elements.
        _units : bool
            Set True if elements are quantities, i.e. are carrying units (default behavior).
            Set to False if elements are raw values, i.e. not carrying units.
        """

        self._attractor = attractor
        self._elements = (
            self._import_elements(*elements) if _units else elements
        )
        self._plane = plane

    @abstractmethod
    def _import_elements(self, *elements):
        """Strip units from elements"""

    @property
    def plane(self):
        """Fundamental plane of the frame."""
        return self._plane

    @property
    def attractor(self):
        """Main attractor."""
        return self._attractor

    @cached_property
    def n(self):
        """Mean motion."""
        return mean_motion(self.attractor.k, self.to_classical().a)

    @cached_property
    def period(self):
        """Period of the orbit."""
        return period(self.attractor.k, self.to_classical().a)

    @cached_property
    def r_p(self):
        """Radius of pericenter."""
        return self.to_classical().a * (1 - self.to_classical().ecc)

    @cached_property
    def r_a(self):
        """Radius of apocenter."""
        return self.to_classical().a * (1 + self.to_classical().ecc)

    @cached_property
    def t_p(self):
        """Elapsed time since latest perifocal passage."""
        return t_p(
            self.to_classical().nu,
            self.to_classical().ecc,
            self.attractor.k,
            self.r_p,
        )

    @abstractmethod
    def to_tuple(self):
        """Expose tuple of values with units"""

    def to_value(self):
        """Expose tuple of values without units (raw)"""
        return self._elements

    def to_vectors(self):
        """Converts to position and velocity vector representation.

        Returns
        -------
        RVState

        """
        raise NotImplementedError

    def to_classical(self):
        """Converts to classical orbital elements representation.

        Returns
        -------
        ClassicalState

        """
        raise NotImplementedError

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation.

        Returns
        -------
        ModifiedEquinoctialState

        """
        raise NotImplementedError


class ClassicalState(BaseState):
    """State defined by its classical orbital elements.

    Orbital elements:

    p : ~astropy.units.Quantity
        Semilatus rectum.
    ecc : ~astropy.units.Quantity
        Eccentricity.
    inc : ~astropy.units.Quantity
        Inclination.
    raan : ~astropy.units.Quantity
        Right ascension of the ascending node.
    argp : ~astropy.units.Quantity
        Argument of the perigee.
    nu : ~astropy.units.Quantity
        True anomaly.

    """

    def _import_elements(self, *elements):
        """Strip units from elements"""
        return (
            elements[0].to_value(u.km),  # p
            elements[1].value,  # ecc
            elements[2].to_value(u.rad),  # inc
            elements[3].to_value(u.rad),  # raan
            elements[4].to_value(u.rad),  # argp
            elements[5].to_value(u.rad),  # nu
        )

    @property
    def p(self):
        """Semilatus rectum."""
        return self._elements[0] << u.km

    @property
    def a(self):
        """Semimajor axis."""
        return self._elements[0] / (1 - self._elements[1] ** 2) << u.km

    @property
    def ecc(self):
        """Eccentricity."""
        return self._elements[1] << u.one

    @property
    def inc(self):
        """Inclination."""
        return self._elements[2] << u.rad

    @property
    def raan(self):
        """Right ascension of the ascending node."""
        return self._elements[3] << u.rad

    @property
    def argp(self):
        """Argument of the perigee."""
        return self._elements[4] << u.rad

    @property
    def nu(self):
        """True anomaly."""
        return self._elements[5] << u.rad

    def to_tuple(self):
        """Expose tuple of values with units"""
        return (
            self.p,
            self.ecc,
            self.inc,
            self.raan,
            self.argp,
            self.nu,
        )

    def to_vectors(self):
        """Converts to position and velocity vector representation."""
        r, v = coe2rv(
            self.attractor.k.to_value(u.km**3 / u.s**2), *self._elements
        )

        return RVState(
            self.attractor,
            (r, v),
            self.plane,
            _units=False,
        )

    def to_classical(self):
        """Converts to classical orbital elements representation."""
        return self

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation."""
        p, f, g, h, k, L = coe2mee(*self._elements)

        return ModifiedEquinoctialState(
            self.attractor,
            (
                p,
                f,
                g,
                h,
                k,
                L,
            ),
            self.plane,
            _units=False,
        )


class RVState(BaseState):
    """State defined by its position and velocity vectors.

    Orbital elements:

    r : ~astropy.units.Quantity
        Position vector wrt attractor center.
    v : ~astropy.units.Quantity
        Velocity vector.

    """

    def _import_elements(self, *elements):
        """Strip units from elements"""
        return (
            elements[0].to_value(u.km),  # r
            elements[1].to_value(u.km / u.s),  # v
        )

    @property
    def r(self):
        """Position vector."""
        return self._elements[0] << u.km

    @property
    def v(self):
        """Velocity vector."""
        return self._elements[1] << u.km / u.s

    def to_tuple(self):
        """Expose tuple of values with units"""
        return (
            self.r,
            self.v,
        )

    def to_vectors(self):
        """Converts to position and velocity vector representation."""
        return self

    def to_classical(self):
        """Converts to classical orbital elements representation."""
        (p, ecc, inc, raan, argp, nu) = rv2coe(
            self.attractor.k.to_value(u.km**3 / u.s**2),
            *self._elements,
        )

        return ClassicalState(
            self.attractor,
            (
                p,
                ecc,
                inc,
                raan,
                argp,
                nu,
            ),
            self.plane,
            _units=False,
        )


class ModifiedEquinoctialState(BaseState):
    """State defined by modified equinoctial elements representation.

    Orbital elements:

    p : ~astropy.units.Quantity
        Semilatus rectum.
    f : ~astropy.units.Quantity
        Second modified equinoctial element.
    g : ~astropy.units.Quantity
        Third modified equinoctial element.
    h : ~astropy.units.Quantity
        Fourth modified equinoctial element.
    k : ~astropy.units.Quantity
        Fifth modified equinoctial element.
    L : ~astropy.units.Quantity
        True longitude.

    """

    def _import_elements(self, *elements):
        """Strip units from elements"""
        return (
            elements[0].to_value(u.km),  # p
            elements[1].to_value(u.rad),  # f
            elements[2].to_value(u.rad),  # g
            elements[3].to_value(u.rad),  # h
            elements[4].to_value(u.rad),  # k
            elements[5].to_value(u.rad),  # L
        )

    @property
    def p(self):
        """Semilatus rectum."""
        return self._elements[0] << u.km

    @property
    def f(self):
        """Second modified equinoctial element."""
        return self._elements[1] << u.rad

    @property
    def g(self):
        """Third modified equinoctial element."""
        return self._elements[2] << u.rad

    @property
    def h(self):
        """Fourth modified equinoctial element."""
        return self._elements[3] << u.rad

    @property
    def k(self):
        """Fifth modified equinoctial element."""
        return self._elements[4] << u.rad

    @property
    def L(self):
        """True longitude."""
        return self._elements[5] << u.rad

    def to_tuple(self):
        """Expose tuple of values with units"""
        return (
            self.p,
            self.f,
            self.g,
            self.h,
            self.k,
            self.L,
        )

    def to_classical(self):
        """Converts to classical orbital elements representation."""
        p, ecc, inc, raan, argp, nu = mee2coe(*self._elements)

        return ClassicalState(
            self.attractor,
            (
                p,
                ecc,
                inc,
                raan,
                argp,
                nu,
            ),
            self.plane,
            _units=False,
        )

    def to_vectors(self):
        """Converts to position and velocity vector representation."""
        r, v = mee2rv(*self._elements)
        return RVState(
            self.attractor,
            (r, v),
            self.plane,
            _units=False,
        )
