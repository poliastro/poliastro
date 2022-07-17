from warnings import warn

from astropy import units as u
from astropy.coordinates import (
    ICRS,
    BarycentricMeanEcliptic,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric_posvel,
)
from astropy.time import Time
from astroquery.jplhorizons import Horizons

from poliastro._math.interpolate import interp1d, sinc_interp, spline_interp
from poliastro.bodies import Earth
from poliastro.frames import Planes
from poliastro.frames.util import get_frame
from poliastro.twobody.sampling import EpochsArray
from poliastro.util import time_range
from poliastro.warnings import TimeScaleWarning

EPHEM_FORMAT = "Ephemerides at {num} epochs from {start} ({start_scale}) to {end} ({end_scale})"


def build_ephem_interpolant(body, period, t_span, rtol=1e-5, attractor=Earth):
    """Interpolates ephemerides data

    Parameters
    ----------
    body : Body
        Source body.
    period : ~astropy.units.Quantity
        Orbital period.
    t_span : list(~astropy.units.Quantity)
        Initial and final epochs.
    rtol : float, optional
        Relative tolerance. Controls the number of sampled data points,
        defaults to 1e-5.
    attractor : ~poliastro.bodies.Body, optional
        Attractor, default to Earth.

    Returns
    -------
    interpolant : callable
        Interpolant function that receives time increment in seconds
        since the initial epoch.

    """
    epochs = time_range(
        Time(t_span[0], format="jd", scale="tdb"),
        end=Time(t_span[1], format="jd", scale="tdb"),
        periods=int(
            ((t_span[1] - t_span[0]) / (period * rtol)).to_value(u.one)
        ),
    )
    ephem = Ephem.from_body(body, epochs, attractor=attractor)

    interpolant = interp1d(
        (epochs - epochs[0]).to_value(u.s),
        ephem._coordinates.xyz.to_value(u.km),
    )
    return interpolant


class BaseInterpolator:
    def interpolate(self, epochs, reference_epochs, coordinates):
        raise NotImplementedError


class SincInterpolator:
    def interpolate(self, epochs, reference_epochs, coordinates):
        def _interp_1d(arr):
            return sinc_interp(arr, reference_epochs.jd, epochs.jd)

        xyz_unit = coordinates.xyz.unit
        d_xyz_unit = coordinates.differentials["s"].d_xyz.unit

        x = _interp_1d(coordinates.x.value) << xyz_unit
        y = _interp_1d(coordinates.y.value) << xyz_unit
        z = _interp_1d(coordinates.z.value) << xyz_unit

        d_x = (
            _interp_1d(coordinates.differentials["s"].d_x.value) << d_xyz_unit
        )
        d_y = (
            _interp_1d(coordinates.differentials["s"].d_y.value) << d_xyz_unit
        )
        d_z = (
            _interp_1d(coordinates.differentials["s"].d_z.value) << d_xyz_unit
        )

        return CartesianRepresentation(
            x, y, z, differentials=CartesianDifferential(d_x, d_y, d_z)
        )


class SplineInterpolator:
    def __init__(self, kind="cubic"):
        self._kind = kind

    def interpolate(self, epochs, reference_epochs, coordinates):
        xyz_unit = coordinates.xyz.unit
        d_xyz_unit = coordinates.differentials["s"].d_xyz.unit

        result_xyz = (
            spline_interp(
                coordinates.xyz.value,
                reference_epochs.jd,
                epochs.jd,
                kind=self._kind,
            )
            << xyz_unit
        )
        result_d_xyz = (
            spline_interp(
                coordinates.differentials["s"].d_xyz.value,
                reference_epochs.jd,
                epochs.jd,
                kind=self._kind,
            )
            << d_xyz_unit
        )

        return CartesianRepresentation(
            result_xyz, differentials=CartesianDifferential(result_d_xyz)
        )


def _get_destination_frame(attractor, plane, epochs):
    if attractor is not None:
        destination_frame = get_frame(attractor, plane, epochs)
    elif plane is Planes.EARTH_ECLIPTIC:
        destination_frame = BarycentricMeanEcliptic()
    else:
        destination_frame = None

    return destination_frame


class Ephem:
    """Time history of position and velocity of some object at particular epochs.

    Instead of creating Ephem objects directly, use the available classmethods.

    Parameters
    ----------
    coordinates : astropy.coordinates.CartesianRepresentation
        Coordinates with velocities.
    epochs : astropy.time.Time
        Epochs corresponding to the coordinates.
    plane : ~poliastro.frames.Planes
        Reference plane of the coordinates.

    """

    def __init__(self, coordinates, epochs, plane):
        if coordinates.ndim != 1 or epochs.ndim != 1:
            raise ValueError(
                f"Coordinates and epochs must have dimension 1, got {coordinates.ndim} and {epochs.ndim}"
            )

        self._epochs = epochs
        self._coordinates = coordinates
        self._plane = Planes(plane)

    def __str__(self):
        return EPHEM_FORMAT.format(
            num=len(self.epochs),
            start=self.epochs[0],
            start_scale=self.epochs[0].scale.upper(),
            end=self.epochs[-1],
            end_scale=self.epochs[-1].scale.upper(),
        )

    def __repr__(self):
        return self.__str__()

    @property
    def epochs(self):
        """Epochs at which the ephemerides was originally sampled."""
        return self._epochs

    @property
    def plane(self):
        """Reference plane of the ephemerides."""
        return self._plane

    @classmethod
    def from_body(
        cls, body, epochs, *, attractor=None, plane=Planes.EARTH_EQUATOR
    ):
        """Return `Ephem` for a `SolarSystemPlanet` at certain epochs.

        Parameters
        ----------
        body : ~poliastro.bodies.SolarSystemPlanet
            Body.
        epochs : ~astropy.time.Time
            Epochs to sample the body positions.
        attractor : ~poliastro.bodies.SolarSystemPlanet, optional
            Body to use as central location,
            if not given the Solar System Barycenter will be used.
        plane : ~poliastro.frames.Planes, optional
            Fundamental plane of the frame, default to Earth Equator.

        """
        if epochs.isscalar:
            epochs = epochs.reshape(1)

        if epochs.scale != "tdb":
            epochs = epochs.tdb
            warn(
                "Input time was converted to scale='tdb' with value "
                f"{epochs.tdb.value}. Use Time(..., scale='tdb') instead.",
                TimeScaleWarning,
                stacklevel=2,
            )

        r, v = get_body_barycentric_posvel(body.name, epochs)
        coordinates = r.with_differentials(
            v.represent_as(CartesianDifferential)
        )

        destination_frame = _get_destination_frame(attractor, plane, epochs)

        if destination_frame is not None:
            coordinates = (
                ICRS(coordinates)
                .transform_to(destination_frame)
                .represent_as(CartesianRepresentation, CartesianDifferential)
            )

        return cls(coordinates, epochs, plane)

    @classmethod
    def from_horizons(
        cls,
        name,
        epochs,
        *,
        attractor=None,
        plane=Planes.EARTH_EQUATOR,
        id_type=None,
    ):
        """Return `Ephem` for an object using JPLHorizons module of Astroquery.

        Parameters
        ----------
        name : str
            Name of the body to query for.
        epochs : ~astropy.time.Time
            Epochs to sample the body positions.
        attractor : ~poliastro.bodies.SolarSystemPlanet, optional
            Body to use as central location,
            if not given the Solar System Barycenter will be used.
        plane : ~poliastro.frames.Planes, optional
            Fundamental plane of the frame, default to Earth Equator.
        id_type : NoneType or str, optional
            Use "smallbody" for Asteroids and Comets and None (default) to first
            search for Planets and Satellites.

        """
        if epochs.isscalar:
            epochs = epochs.reshape(1)

        refplanes_dict = {
            Planes.EARTH_EQUATOR: "earth",
            Planes.EARTH_ECLIPTIC: "ecliptic",
        }
        refplane = refplanes_dict[plane]

        if attractor is not None:
            bodies_dict = {
                "sun": 10,
                "mercury": 199,
                "venus": 299,
                "earth": 399,
                "mars": 499,
                "jupiter": 599,
                "saturn": 699,
                "uranus": 799,
                "neptune": 899,
            }
            location = f"500@{bodies_dict[attractor.name.lower()]}"
        else:
            location = "@ssb"

        obj = Horizons(
            id=name, location=location, epochs=epochs.jd, id_type=id_type
        ).vectors(refplane=refplane)

        x = obj["x"]
        y = obj["y"]
        z = obj["z"]
        d_x = obj["vx"]
        d_y = obj["vy"]
        d_z = obj["vz"]

        coordinates = CartesianRepresentation(
            x, y, z, differentials=CartesianDifferential(d_x, d_y, d_z)
        )
        return cls(coordinates, epochs, plane)

    @classmethod
    def from_orbit(
        cls,
        orbit,
        epochs,
        plane=Planes.EARTH_EQUATOR,
    ):
        """Return `Ephem` from an `Orbit` at certain epochs.

        Parameters
        ----------
        orbit : ~poliastro.twobody.orbit.Orbit
            Orbit.
        epochs : ~astropy.time.Time
            Epochs to sample the orbit positions.
        plane : ~poliastro.frames.Planes, optional
            Fundamental plane of the frame, default to Earth Equator.

        """
        if epochs.isscalar:
            epochs = epochs.reshape(1)

        return orbit.change_plane(plane).to_ephem(strategy=EpochsArray(epochs))

    def sample(self, epochs=None, *, interpolator=SplineInterpolator()):
        """Returns coordinates at specified epochs.

        Parameters
        ----------
        epochs : ~astropy.time.Time, optional
            Epochs to sample the ephemerides,
            if not given the original one from the object will be used.
        interpolator : ~poliastro.ephem.BaseInterpolator, optional
            Interpolation method to use for epochs outside of the original ones,
            default to splines.

        Returns
        -------
        CartesianRepresentation
            Sampled coordinates with velocities.

        """
        if epochs is None or epochs.isscalar and (epochs == self.epochs).all():
            return self._coordinates

        coordinates = interpolator.interpolate(
            epochs.reshape(-1),
            self.epochs,
            self._coordinates,
        )

        return coordinates

    def rv(self, epochs=None, **kwargs):
        """Position and velocity vectors at given epochs.

        Parameters
        ----------
        epochs : ~astropy.time.Time, optional
            Epochs to sample the ephemerides, default to now.
        **kwargs
            Extra kwargs for interpolation method.

        """
        coordinates = self.sample(epochs, **kwargs)

        r = coordinates.get_xyz(xyz_axis=1)
        v = coordinates.differentials["s"].get_d_xyz(xyz_axis=1)

        if epochs is not None and epochs.isscalar:
            return r[0], v[0]
        else:
            return r, v
