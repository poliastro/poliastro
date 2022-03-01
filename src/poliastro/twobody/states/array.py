from abc import ABC, abstractmethod

from astropy import time, units as u
import numpy as np

from poliastro.core.elements import coe2mee, coe2rv, mee2coe, mee2rv, rv2coe
from poliastro.twobody.elements import mean_motion, period

from .scalar import (
    ClassicalState,
    ModifiedEquinoctialState,
    RVState,
)


class BaseStateArray(ABC):
    """Base State class, meant to be subclassed."""

    @abstractmethod
    def __init__(self, epoch, attractor, plane):
        """Constructor.

        Parameters
        ----------
        epoch : ~astropy.time.Time
            Array of epochs of each individual state (orbit).
        attractor : Body
            Common main attractor.
        plane : ~poliastro.frames.enums.Planes
            Common reference plane for the elements.

        """
        self._epoch = epoch  # type: time.Time
        self._attractor = attractor
        self._plane = plane

    @abstractmethod
    def __getitem__(self, idx):
        """Get item or slice from state array."""
        raise NotImplementedError

    @abstractmethod
    def __setitem__(self, idx, value):
        """Set item or slice from state array."""
        raise NotImplementedError

    @abstractmethod
    def copy(self):
        """Copy state array."""
        raise NotImplementedError

    @abstractmethod
    def reshape(self, *args):
        """Reshape state array."""
        raise NotImplementedError

    @property
    def ndim(self):
        """Number of dimensions of array."""
        return self._epoch.ndim

    @property
    def shape(self):
        """Shape of array."""
        return self._epoch.shape

    @property
    def size(self):
        """Size of array."""
        return self._epoch.size

    @property
    def epoch(self):
        """Array of epochs of each individual state (orbit)."""
        return self._epoch

    @property
    def plane(self):
        """Common fundamental plane of the frame."""
        return self._plane

    @property
    def attractor(self):
        """Common main attractor."""
        return self._attractor

    @property
    def n(self):
        """Mean motion array."""
        return mean_motion(self.attractor.k, self.to_classical().a)

    @property
    def period(self):
        """Period of the orbit array."""
        return period(self.attractor.k, self.to_classical().a)

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


class ClassicalStateArray(BaseStateArray):
    """State defined by its classical orbital elements."""

    def __init__(self, epoch, attractor, p, ecc, inc, raan, argp, nu, plane):
        """Constructor.

        Parameters
        ----------
        epoch : ~astropy.time.Time
            Array of epochs of each individual state (orbit).
        attractor : Body
            Common main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum array.
        ecc : ~astropy.units.Quantity
            Eccentricity array.
        inc : ~astropy.units.Quantity
            Inclination array.
        raan : ~astropy.units.Quantity
            Right ascension of the ascending node array.
        argp : ~astropy.units.Quantity
            Argument of the perigee array.
        nu : ~astropy.units.Quantity
            True anomaly array.
        plane : ~poliastro.frames.enums.Planes
            Common reference plane for the elements.
        """
        assert epoch.spahe == p.shape == ecc.shape == inc.shape == raan.shape == argp.shape == nu.shape
        super().__init__(epoch, attractor, plane)
        self._p = p
        self._ecc = ecc
        self._inc = inc
        self._raan = raan
        self._argp = argp
        self._nu = nu

    def __getitem__(self, idx):
        """Get item or slice from state array."""
        epoch = self._epoch[idx]
        p = self._p[idx]
        ecc = self._ecc[idx]
        inc = self._inc[idx]
        raan = self._raan[idx]
        argp = self._argp[idx]
        nu = self._nu[idx]
        cls = ClassicalState if p.ndim == 0 else type(self)
        return cls(
            epoch = epoch,
            attractor = self._attractor,
            p = p,
            ecc = ecc,
            inc = inc,
            raan = raan,
            argp = argp,
            nu = nu,
            plane = self._plane,
        )  # type: ignore

    def __setitem__(self, idx, value):
        """Set item or slice from state array."""
        raise NotImplementedError # TODO

    def copy(self):
        """Copy state array."""
        return type(self)(
            epoch = self._epoch.copy(),
            attractor = self._attractor,
            p = self._p.copy(),
            ecc = self._ecc.copy(),
            inc = self._inc.copy(),
            raan = self._raan.copy(),
            argp = self._argp.copy(),
            nu = self._nu.copy(),
            plane = self._plane,
        )

    def reshape(self, *args):
        """Reshape state array."""
        return type(self)(
            epoch = self._epoch.reshape(*args),
            attractor = self._attractor,
            p = self._p.reshape(*args),
            ecc = self._ecc.reshape(*args),
            inc = self._inc.reshape(*args),
            raan = self._raan.reshape(*args),
            argp = self._argp.reshape(*args),
            nu = self._nu.reshape(*args),
            plane = self._plane,
        )

    @property
    def p(self):
        """Semilatus rectum array."""
        return self._p

    @property
    def a(self):
        """Semimajor axis array."""
        return self.p / (1 - self.ecc**2)

    @property
    def ecc(self):
        """Eccentricity array."""
        return self._ecc

    @property
    def inc(self):
        """Inclination array."""
        return self._inc

    @property
    def raan(self):
        """Right ascension of the ascending node array."""
        return self._raan

    @property
    def argp(self):
        """Argument of the perigee array."""
        return self._argp

    @property
    def nu(self):
        """True anomaly array."""
        return self._nu

    def to_vectors(self):
        """Converts to position and velocity vector representation."""

        r = np.zeros(self._p.shape, dtype = self._p.dtype)
        v = np.zeros(self._p.shape, dtype = self._p.dtype)

        r_flat = r.reshape((-1,))
        v_flat = v.reshape((-1,))

        k = self._attractor.k.to_value(u.km**3 / u.s**2)
        p_flat = self.p.to_value(u.km).reshape((-1,))
        ecc_flat = self.ecc.value.reshape((-1,))
        inc_flat = self.inc.to_value(u.rad).reshape((-1,))
        raan_flat = self.raan.to_value(u.rad).reshape((-1,))
        argp_flat = self.argp.to_value(u.rad).reshape((-1,))
        nu_flat = self.nu.to_value(u.rad).reshape((-1,))

        for idx in range(r_flat.size):  # TODO replace with vector implementation
            r_flat[idx], v_flat[idx] = coe2rv(
                k,
                p_flat[idx],
                ecc_flat[idx],
                inc_flat[idx],
                raan_flat[idx],
                argp_flat[idx],
                nu_flat[idx],
            )

        return RVStateArray(self.epoch.copy(), self.attractor, r * u.km, v * u.km / u.s, self.plane)

    def to_classical(self):
        """Converts to classical orbital elements representation."""
        return self.copy()

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation."""

        p = np.zeros(self._p.shape, dtype = self._p.dtype)
        f = np.zeros(self._p.shape, dtype = self._p.dtype)
        g = np.zeros(self._p.shape, dtype = self._p.dtype)
        h = np.zeros(self._p.shape, dtype = self._p.dtype)
        k = np.zeros(self._p.shape, dtype = self._p.dtype)
        L = np.zeros(self._p.shape, dtype = self._p.dtype)

        p_flat = p.reshape((-1,))
        f_flat = f.reshape((-1,))
        g_flat = g.reshape((-1,))
        h_flat = h.reshape((-1,))
        k_flat = k.reshape((-1,))
        L_flat = L.reshape((-1,))

        p__flat = self.p.to_value(u.km).reshape((-1,))
        ecc_flat = self.ecc.value.reshape((-1,))
        inc_flat = self.inc.to_value(u.rad).reshape((-1,))
        raan_flat = self.raan.to_value(u.rad).reshape((-1,))
        argp_flat = self.argp.to_value(u.rad).reshape((-1,))
        nu_flat = self.nu.to_value(u.rad).reshape((-1,))

        for idx in range(p_flat.size):  # TODO replace with vector implementation
            p_flat[idx], f_flat[idx], g_flat[idx], h_flat[idx], k_flat[idx], L_flat[idx] = coe2mee(
                p__flat[idx],
                ecc_flat[idx],
                inc_flat[idx],
                raan_flat[idx],
                argp_flat[idx],
                nu_flat[idx],
            )

        return ModifiedEquinoctialStateArray(
            self.epoch.copy(),
            self.attractor,
            p * u.km,
            f * u.rad,
            g * u.rad,
            h * u.rad,
            k * u.rad,
            L * u.rad,
            self.plane,
        )


class RVStateArray(BaseStateArray):
    """State defined by its position and velocity vectors."""

    def __init__(self, epoch, attractor, r, v, plane):
        """Constructor.

        Parameters
        ----------
        epoch : ~astropy.time.Time
            Array of epochs of each individual state (orbit).
        attractor : Body
            Common main attractor.
        r : ~astropy.units.Quantity
            Position vector wrt attractor center array.
        v : ~astropy.units.Quantity
            Velocity vector array.
        plane : ~poliastro.frames.enums.Planes
            Common reference plane for the elements.

        """
        assert epoch.shape == r.shape == v.shape
        super().__init__(epoch, attractor, plane)
        self._r = r
        self._v = v

    def __getitem__(self, idx):
        """Get item or slice from state array."""
        epoch = self._epoch[idx]
        r = self._r[idx]
        v = self._v[idx]
        cls = RVState if r.ndim == 0 else type(self)
        return cls(
            epoch = epoch,
            attractor = self._attractor,
            r = r,
            v = v,
            plane = self._plane,
        )  # type: ignore

    def __setitem__(self, idx, value):
        """Set item or slice from state array."""
        raise NotImplementedError # TODO

    def copy(self):
        """Copy state array."""
        return type(self)(
            epoch = self._epoch.copy(),
            attractor = self._attractor,
            r = self._r.copy(),
            v = self._v.copy(),
            plane = self._plane,
        )

    def reshape(self, *args):
        """Reshape state array."""
        return type(self)(
            epoch = self._epoch.reshape(*args),
            attractor = self._attractor,
            r = self._r.reshape(*args),
            v = self._v.reshape(*args),
            plane = self._plane,
        )

    @property
    def r(self):
        """Position vector array."""
        return self._r

    @property
    def v(self):
        """Velocity vector array."""
        return self._v

    def to_vectors(self):
        """Converts to position and velocity vector representation."""
        return self.copy()

    def to_classical(self):
        """Converts to classical orbital elements representation."""

        p = np.zeros(self._r.shape, dtype = self._r.dtype)
        ecc = np.zeros(self._r.shape, dtype = self._r.dtype)
        inc = np.zeros(self._r.shape, dtype = self._r.dtype)
        raan = np.zeros(self._r.shape, dtype = self._r.dtype)
        argp = np.zeros(self._r.shape, dtype = self._r.dtype)
        nu = np.zeros(self._r.shape, dtype = self._r.dtype)

        p_flat = p.reshape((-1,))
        ecc_flat = ecc.reshape((-1,))
        inc_flat = inc.reshape((-1,))
        raan_flat = raan.reshape((-1,))
        argp_flat = argp.reshape((-1,))
        nu_flat = nu.reshape((-1,))

        k = self.attractor.k.to_value(u.km**3 / u.s**2)
        r_flat = self.r.to_value(u.km).reshape((-1,))
        v_flat = self.v.to_value(u.km / u.s).reshape((-1,))

        for idx in range(p_flat.size):  # TODO replace with vector implementation
            p_flat[idx], ecc_flat[idx], inc_flat[idx], raan_flat[idx], argp_flat[idx], nu_flat[idx] = rv2coe(
                k,
                r_flat[idx],
                v_flat[idx],
            )

        return ClassicalStateArray(
            self.epoch.copy(),
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
            self.plane,
        )


class ModifiedEquinoctialStateArray(BaseStateArray):
    """State defined by modified equinoctial elements representation."""

    def __init__(self, epoch, attractor, p, f, g, h, k, L, plane):
        """Constructor.

        Parameters
        ----------
        epoch : ~astropy.time.Time
            Array of epochs of each individual state (orbit).
        attractor : Body
            Common main attractor.
        p : ~astropy.units.Quantity
            Semilatus rectum array.
        f : ~astropy.units.Quantity
            Second modified equinoctial element array.
        g : ~astropy.units.Quantity
            Third modified equinoctial element array.
        h : ~astropy.units.Quantity
            Fourth modified equinoctial element array.
        k : ~astropy.units.Quantity
            Fifth modified equinoctial element array.
        L : ~astropy.units.Quantity
            True longitude array.
        plane : ~poliastro.frames.enums.Planes
            Common reference plane for the elements.

        """
        assert epoch.shape == p.shape == f.shape == g.shape == h.shape == k.shape == L.shape
        super().__init__(epoch, attractor, plane)
        self._p = p
        self._f = f
        self._g = g
        self._h = h
        self._k = k
        self._L = L

    def __getitem__(self, idx):
        """Get item or slice from state array."""
        epoch = self._epoch[idx]
        p = self._p[idx]
        f = self._f[idx]
        g = self._g[idx]
        h = self._h[idx]
        k = self._k[idx]
        L = self._L[idx]
        cls = ModifiedEquinoctialState if p.ndim == 0 else type(self)
        return cls(
            epoch = epoch,
            attractor = self._attractor,
            p = p,
            f = f,
            g = g,
            h = h,
            k = k,
            L = L,
            plane = self._plane,
        )  # type: ignore

    def __setitem__(self, idx, value):
        """Set item or slice from state array."""
        raise NotImplementedError # TODO

    def copy(self):
        """Copy state array."""
        return type(self)(
            epoch = self._epoch.copy(),
            attractor = self._attractor,
            p = self._p.copy(),
            f = self._f.copy(),
            g = self._g.copy(),
            h = self._h.copy(),
            k = self._k.copy(),
            L = self._L.copy(),
            plane = self._plane,
        )

    def reshape(self, *args):
        """Reshape state array."""
        return type(self)(
            epoch = self._epoch.reshape(*args),
            attractor = self._attractor,
            p = self._p.reshape(*args),
            f = self._f.reshape(*args),
            g = self._g.reshape(*args),
            h = self._h.reshape(*args),
            k = self._k.reshape(*args),
            L = self._L.reshape(*args),
            plane = self._plane,
        )

    @property
    def p(self):
        """Semilatus rectum array."""
        return self._p

    @property
    def f(self):
        """Second modified equinoctial element array."""
        return self._f

    @property
    def g(self):
        """Third modified equinoctial element array."""
        return self._g

    @property
    def h(self):
        """Fourth modified equinoctial element array."""
        return self._h

    @property
    def k(self):
        """Fifth modified equinoctial element array."""
        return self._k

    @property
    def L(self):
        """True longitude array."""
        return self._L

    def to_classical(self):
        """Converts to classical orbital elements representation."""

        p = np.zeros(self._p.shape, dtype = self._p.dtype)
        ecc = np.zeros(self._p.shape, dtype = self._p.dtype)
        inc = np.zeros(self._p.shape, dtype = self._p.dtype)
        raan = np.zeros(self._p.shape, dtype = self._p.dtype)
        argp = np.zeros(self._p.shape, dtype = self._p.dtype)
        nu = np.zeros(self._p.shape, dtype = self._p.dtype)

        p_flat = p.reshape((-1,))
        ecc_flat = ecc.reshape((-1,))
        inc_flat = inc.reshape((-1,))
        raan_flat = raan.reshape((-1,))
        argp_flat = argp.reshape((-1,))
        nu_flat = nu.reshape((-1,))

        p__flat = self.p.to_value(u.km).reshape((-1,))
        f_flat = self.f.to_value(u.rad).reshape((-1,))
        g_flat = self.g.to_value(u.rad).reshape((-1,))
        h_flat = self.h.to_value(u.rad).reshape((-1,))
        k_flat = self.k.to_value(u.rad).reshape((-1,))
        L_flat = self.L.to_value(u.rad).reshape((-1,))

        for idx in range(p_flat.size):  # TODO replace with vector implementation
            p_flat[idx], ecc_flat[idx], inc_flat[idx], raan_flat[idx], argp_flat[idx], nu_flat[idx] = mee2coe(
                p__flat[idx],
                f_flat[idx],
                g_flat[idx],
                h_flat[idx],
                k_flat[idx],
                L_flat[idx],
            )

        return ClassicalStateArray(
            self.epoch.copy(),
            self.attractor,
            p * u.km,
            ecc * u.one,
            inc * u.rad,
            raan * u.rad,
            argp * u.rad,
            nu * u.rad,
            self.plane,
        )

    def to_vectors(self):
        """Converts to position and velocity vector representation."""

        r = np.zeros(self._p.shape, dtype = self._p.dtype)
        v = np.zeros(self._p.shape, dtype = self._p.dtype)

        r_flat = r.reshape((-1,))
        v_flat = v.reshape((-1,))

        p_flat = self.p.to(u.km).value.reshape((-1,))
        f_flat = self.f.to(u.rad).value.reshape((-1,))
        g_flat = self.g.to(u.rad).value.reshape((-1,))
        h_flat = self.h.to(u.rad).value.reshape((-1,))
        k_flat = self.k.to(u.rad).value.reshape((-1,))
        L_flat = self.L.to(u.rad).value.reshape((-1,))

        for idx in range(r_flat.size):  # TODO replace with vector implementation
            r_flat[idx], v_flat[idx] = mee2rv(
                p_flat[idx],
                f_flat[idx],
                g_flat[idx],
                h_flat[idx],
                k_flat[idx],
                L_flat[idx],
            )

        return RVStateArray(self.epoch.copy(), self.attractor, r * u.km, v * u.km / u.s, self.plane)

    def to_equinoctial(self):
        """Converts to modified equinoctial elements representation."""
        return self.copy()
