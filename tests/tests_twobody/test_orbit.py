import pickle

import matplotlib
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import (
    ITRS,
    CartesianDifferential,
    CartesianRepresentation,
    SkyCoord,
)
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from numpy.testing import assert_allclose, assert_array_equal

from poliastro.bodies import (
    Body,
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
)
from poliastro.constants import J2000, J2000_TDB
from poliastro.examples import iss
from poliastro.frames.ecliptic import HeliocentricEclipticJ2000
from poliastro.frames.enums import Planes
from poliastro.frames.equatorial import (
    GCRS,
    HCRS,
    ICRS,
    JupiterICRS,
    MarsICRS,
    MercuryICRS,
    NeptuneICRS,
    PlutoICRS,
    SaturnICRS,
    UranusICRS,
    VenusICRS,
)
from poliastro.frames.util import get_frame
from poliastro.twobody.orbit import (
    Orbit,
    OrbitSamplingWarning,
    PatchedConicsWarning,
    TimeScaleWarning,
)


@pytest.fixture()
def hyperbolic():
    r = [1.197659243752796e09, -4.443716685978071e09, -1.747610548576734e09] * u.km
    v = (
        [5.540549267188614e00, -1.251544669134140e01, -4.848892572767733e00]
        * u.km
        / u.s
    )
    epoch = Time("2015-07-14 07:59", scale="tdb")
    return Orbit.from_vectors(Sun, r, v, epoch)


def test_default_time_for_new_state():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    expected_epoch = J2000
    ss = Orbit.from_classical(_body, _d, _, _a, _a, _a, _a)
    assert ss.epoch == expected_epoch


def test_state_raises_unitserror_if_elements_units_are_wrong():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    wrong_angle = 1.0 * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        Orbit.from_classical(Sun, _d, _, _a, _a, _a, wrong_angle)
    assert (
        "UnitsError: Argument 'nu' to function 'from_classical' must be in units convertible to 'rad'."
        in excinfo.exconly()
    )


def test_state_raises_unitserror_if_rv_units_are_wrong():
    _d = [1.0, 0.0, 0.0] * u.AU
    wrong_v = [0.0, 1.0e-6, 0.0] * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        Orbit.from_vectors(Sun, _d, wrong_v)
    assert (
        "UnitsError: Argument 'v' to function 'from_vectors' must be in units convertible to 'm / s'."
        in excinfo.exconly()
    )


def test_parabolic_elements_fail_early():
    attractor = Earth
    ecc = 1.0 * u.one
    _d = 1.0 * u.AU  # Unused distance
    _a = 1.0 * u.deg  # Unused angle
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_classical(attractor, _d, ecc, _a, _a, _a, _a)
    assert (
        "ValueError: For parabolic orbits use Orbit.parabolic instead"
        in excinfo.exconly()
    )


def test_bad_inclination_raises_exception():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    bad_inc = 200 * u.deg
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_classical(_body, _d, _, bad_inc, _a, _a, _a)
    assert (
        "ValueError: Inclination must be between 0 and 180 degrees" in excinfo.exconly()
    )


def test_bad_hyperbolic_raises_exception():
    bad_a = 1.0 * u.AU
    ecc = 1.5 * u.one
    _inc = 100 * u.deg  # Unused inclination
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_classical(_body, bad_a, ecc, _inc, _a, _a, _a)
    assert "Hyperbolic orbits have negative semimajor axis" in excinfo.exconly()


def test_apply_maneuver_changes_epoch():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = Orbit.from_classical(Sun, _d, _, _a, _a, _a, _a)
    dt = 1 * u.h
    dv = [0, 0, 0] * u.km / u.s
    orbit_new = ss.apply_maneuver([(dt, dv)])
    assert orbit_new.epoch == ss.epoch + dt


def test_orbit_from_ephem_with_no_epoch_is_today():
    # This is not that obvious http://stackoverflow.com/q/6407362/554319
    body = Earth
    ss = Orbit.from_body_ephem(body)
    assert (Time.now() - ss.epoch).sec < 1


def test_from_ephem_raises_warning_if_time_is_not_tdb_with_proper_time(recwarn):
    body = Earth
    epoch = Time("2017-09-29 07:31:26", scale="utc")
    expected_epoch_string = "2017-09-29 07:32:35.182"  # epoch.tdb.value

    Orbit.from_body_ephem(body, epoch)

    w = recwarn.pop(TimeScaleWarning)
    assert expected_epoch_string in str(w.message)


@pytest.mark.parametrize("body", [Moon, Pluto])
def test_from_ephem_raises_error_for_pluto_moon(body):
    with pytest.raises(RuntimeError) as excinfo:
        Orbit.from_body_ephem(body)
    assert "To compute the position and velocity" in excinfo.exconly()


def test_circular_has_proper_semimajor_axis():
    alt = 500 * u.km
    attractor = Earth
    expected_a = Earth.R + alt
    ss = Orbit.circular(attractor, alt)
    assert ss.a == expected_a


def test_geosync_has_proper_period():
    expected_period = 1436 * u.min

    ss = Orbit.circular(Earth, alt=42164 * u.km - Earth.R)

    assert_quantity_allclose(ss.period, expected_period, rtol=1e-4)


def test_parabolic_has_proper_eccentricity():
    attractor = Earth
    _d = 1.0 * u.AU  # Unused distance
    _a = 1.0 * u.deg  # Unused angle
    expected_ecc = 1.0 * u.one
    ss = Orbit.parabolic(attractor, _d, _a, _a, _a, _a)
    assert_allclose(ss.ecc, expected_ecc)


def test_parabolic_has_zero_energy():
    attractor = Earth
    _d = 1.0 * u.AU  # Unused distance
    _a = 1.0 * u.deg  # Unused angle
    ss = Orbit.parabolic(attractor, _d, _a, _a, _a, _a)
    assert_allclose(ss.energy.value, 0.0, atol=1e-16)


def test_pqw_for_circular_equatorial_orbit():
    ss = Orbit.circular(Earth, 600 * u.km)
    expected_p = [1, 0, 0] * u.one
    expected_q = [0, 1, 0] * u.one
    expected_w = [0, 0, 1] * u.one
    p, q, w = ss.pqw()
    assert_allclose(p, expected_p)
    assert_allclose(q, expected_q)
    assert_allclose(w, expected_w)


@pytest.mark.parametrize(
    "attractor,alt,argp,expected_argp,expected_inc",
    [
        (
            Earth,
            1e6 * u.m,
            3 * np.pi / 2 * u.rad,
            3 * np.pi / 2 * u.rad,
            63.4349 * np.pi / 180 * u.rad,
        ),
        (Mars, 3e8 * u.m, 0 * u.deg, 0 * u.deg, 63.4349 * np.pi / 180 * u.rad),
    ],
)
def test_frozen_orbit_argp(attractor, alt, argp, expected_argp, expected_inc):
    orbit = Orbit.frozen(attractor, alt, argp=argp)
    assert_allclose(orbit.argp, expected_argp)
    assert_allclose(orbit.inc, expected_inc)


@pytest.mark.parametrize(
    "attractor,alt,inc,argp,expected_inc,expected_argp",
    [
        (
            Earth,
            1e6 * u.m,
            116.5651 * np.pi / 180 * u.rad,
            3 * np.pi / 2 * u.rad,
            116.5651 * np.pi / 180 * u.rad,
            3 * np.pi / 2 * u.rad,
        ),
        (
            Mars,
            3e8 * u.m,
            63.4349 * np.pi / 180 * u.rad,
            np.pi / 2 * u.rad,
            63.4349 * np.pi / 180 * u.rad,
            np.pi / 2 * u.rad,
        ),
    ],
)
def test_frozen_orbit_with_critical_argp_and_critical_inc(
    attractor, alt, inc, argp, expected_inc, expected_argp
):
    orbit = Orbit.frozen(attractor, alt, inc=inc, argp=argp)
    assert_allclose(orbit.argp, expected_argp)
    assert_allclose(orbit.inc, expected_inc)


@pytest.mark.parametrize(
    "attractor,alt,expected_inc,expected_argp",
    [
        (Earth, 1e6 * u.m, 63.4349 * np.pi / 180 * u.rad, np.pi / 2 * u.rad),
        (Mars, 3e8 * u.m, 63.4349 * np.pi / 180 * u.rad, np.pi / 2 * u.rad),
    ],
)
def test_frozen_orbit_no_args(attractor, alt, expected_inc, expected_argp):
    orbit = Orbit.frozen(attractor, alt)
    argp = orbit.argp
    inc = orbit.inc
    assert_allclose(argp, expected_argp)
    assert_allclose(inc, expected_inc)


@pytest.mark.parametrize(
    "attractor,alt,argp,expected_inc,ecc,expected_ecc",
    [
        (
            Earth,
            1e6 * u.m,
            2 * u.deg,  # Non critical value
            63.4349 * np.pi / 180 * u.rad,
            None,
            0.0549 * u.one,
        ),
        (
            Mars,
            3e8 * u.m,
            0 * u.deg,  # Non critical value
            63.4349 * np.pi / 180 * u.rad,
            0.04 * u.one,
            0.04 * u.one,
        ),
    ],
)
def test_frozen_orbit_with_non_critical_argp(
    attractor, alt, argp, expected_inc, ecc, expected_ecc
):
    orbit = Orbit.frozen(attractor, alt, argp=argp, ecc=ecc)  # Non-critical value
    assert_allclose(orbit.inc, expected_inc)
    assert_allclose(orbit.ecc, expected_ecc)


def test_frozen_orbit_non_critical_inclination():
    orbit = Orbit.frozen(Earth, 1e3 * u.km, inc=0 * u.deg)  # Non-critical value
    assert orbit.argp in [np.pi / 2, 3 * np.pi / 2] * u.rad


def test_frozen_orbit_venus_special_case():
    with pytest.raises(NotImplementedError) as excinfo:
        Orbit.frozen(Venus, 1 * u.m)
    assert excinfo.type == NotImplementedError
    assert str(excinfo.value) == "This has not been implemented for Venus"


def test_frozen_orbit_non_spherical_arguments():
    with pytest.raises(AttributeError) as excinfo:
        Orbit.frozen(Jupiter, 1 * u.m)
    assert excinfo.type == AttributeError
    assert (
        str(excinfo.value)
        == "Attractor Jupiter has not spherical harmonics implemented"
    )


def test_frozen_orbit_altitude():
    with pytest.raises(ValueError) as excinfo:
        Orbit.frozen(Earth, -1 * u.m)
    assert excinfo.type == ValueError
    assert (
        str(excinfo.value)
        == "The semimajor axis may not be smaller that Earth's radius"
    )


def test_orbit_representation():
    ss = Orbit.circular(
        Earth, 600 * u.km, 20 * u.deg, epoch=Time("2018-09-08 09:04:00", scale="tdb")
    )
    expected_str = "6978 x 6978 km x 20.0 deg (GCRS) orbit around Earth (\u2641) at epoch 2018-09-08 09:04:00.000 (TDB)"

    assert str(ss) == repr(ss) == expected_str


def test_orbit_no_frame_representation():
    date_launch = Time("2011-11-26 15:02", scale="utc")
    r = [61_445.76498656, 24_827.93010168, 0.0] * u.km
    v = [-0.42581645, -0.18867869, 0.0] * u.km / u.s
    ss = Orbit.from_vectors(Moon, r, v, date_launch)
    expected_str = "106 x -142299 km x 180.0 deg orbit around Moon (\u263E) at epoch 2011-11-26 15:02:00.000 (UTC)"

    assert str(ss) == repr(ss) == expected_str


def test_sample_numpoints():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    ss = Orbit.from_classical(_body, _d, _, _a, _a, _a, _a)
    positions = ss.sample(values=50)
    assert len(positions) == 50


@pytest.mark.parametrize("num_points", [3, 5, 7, 9, 11, 101])
def test_sample_num_points(num_points):
    # Data from Vallado, example 2.4
    r0 = [1_131.340, -2_282.343, 6_672.423] * u.km
    v0 = [-5.64305, 4.30333, 2.42879] * u.km / u.s
    ss0 = Orbit.from_vectors(Earth, r0, v0)

    # TODO: Test against the perigee and apogee
    # expected_ss = ss0.propagate(ss0.period / 2)

    rr = ss0.sample(num_points)

    assert len(rr) == num_points
    # assert_quantity_allclose(rr[num_points // 2].data.xyz, expected_ss.r)


def test_sample_big_orbits():
    # See https://github.com/poliastro/poliastro/issues/265
    ss = Orbit.from_vectors(
        Sun,
        [-9_018_878.6, -94_116_055, 22_619_059] * u.km,
        [-49.950923, -12.948431, -4.2925158] * u.km / u.s,
    )
    positions = ss.sample(15)
    assert len(positions) == 15


def test_hyperbolic_nu_value_check(hyperbolic):
    positions = hyperbolic.sample(100)

    assert isinstance(positions, CartesianRepresentation)
    assert len(positions) == 100


def test_hyperbolic_modulus_wrapped_nu():
    ss = Orbit.from_vectors(
        Sun,
        [-9.77441841e07, 1.01000539e08, 4.37584668e07] * u.km,
        [23.75936985, -43.09599568, -8.7084724] * u.km / u.s,
    )
    num_values = 3

    positions = ss.sample(num_values)

    assert_quantity_allclose(positions[0].xyz, ss.r)


@pytest.mark.parametrize("min_anomaly", [-30 * u.deg, -10 * u.deg])
@pytest.mark.parametrize("max_anomaly", [10 * u.deg, 30 * u.deg])
def test_sample_hyperbolic_limits(hyperbolic, min_anomaly, max_anomaly):
    num_points = 50

    coords = hyperbolic.sample(
        num_points, min_anomaly=min_anomaly, max_anomaly=max_anomaly
    )

    assert len(coords) == num_points


def test_sample_hyperbolic_outside_limits(hyperbolic):
    with pytest.warns(OrbitSamplingWarning, match="anomaly outside range, clipping"):
        hyperbolic.sample(3, min_anomaly=-np.pi * u.rad)

    with pytest.warns(OrbitSamplingWarning, match="anomaly outside range, clipping"):
        hyperbolic.sample(3, max_anomaly=np.pi * u.rad)


def test_orbit_is_pickable(hyperbolic):
    pickled = pickle.dumps(hyperbolic)
    ss_result = pickle.loads(pickled)

    assert_array_equal(hyperbolic.r, ss_result.r)
    assert_array_equal(hyperbolic.v, ss_result.v)
    assert ss_result.epoch == hyperbolic.epoch


def test_orbit_plot_is_static():
    # Data from Curtis, example 4.3
    r = [-6_045, -3_490, 2_500] * u.km
    v = [-3.457, 6.618, 2.533] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)

    plot = ss.plot()

    assert isinstance(plot[0], matplotlib.lines.Line2D)
    assert isinstance(plot[1], matplotlib.lines.Line2D)


def test_orbit_plot_static_3d():
    # Data from Curtis, example 4.3
    r = [-6_045, -3_490, 2_500] * u.km
    v = [-3.457, 6.618, 2.533] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)
    with pytest.raises(
        ValueError,
        match="The static plotter does not support 3D, use `interactive=True`",
    ):
        ss.plot(use_3d=True)


@pytest.mark.parametrize("use_3d", [False, True])
def test_orbit_plot_is_not_static(use_3d):
    from plotly.graph_objects import Figure

    # Data from Curtis, example 4.3
    r = [-6_045, -3_490, 2_500] * u.km
    v = [-3.457, 6.618, 2.533] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)

    plot = ss.plot(interactive=True, use_3d=use_3d)

    assert isinstance(plot, Figure)


@pytest.mark.parametrize(
    "attractor, expected_frame_class",
    [
        (Sun, HCRS),
        (Mercury, MercuryICRS),
        (Venus, VenusICRS),
        (Earth, GCRS),
        (Mars, MarsICRS),
        (Jupiter, JupiterICRS),
        (Saturn, SaturnICRS),
        (Uranus, UranusICRS),
        (Neptune, NeptuneICRS),
        (Pluto, PlutoICRS),
    ],
)
def test_orbit_get_frame_returns_proper_frame(attractor, expected_frame_class):
    # Dummy data
    r = [1e09, -4e09, -1e09] * u.km
    v = [5e00, -1e01, -4e00] * u.km / u.s
    epoch = Time("2015-07-14 07:59", scale="tdb")

    ss = Orbit.from_vectors(attractor, r, v, epoch)
    frame = ss.get_frame()

    assert frame.is_equivalent_frame(expected_frame_class(obstime=epoch))
    assert frame.obstime == epoch


def test_orbit_from_custom_body_raises_error_when_asked_frame():
    attractor = Body(Sun, 1 * u.km ** 3 / u.s ** 2, "_DummyPlanet")

    r = [1e09, -4e09, -1e09] * u.km
    v = [5e00, -1e01, -4e00] * u.km / u.s

    ss = Orbit.from_vectors(attractor, r, v)

    with pytest.raises(NotImplementedError) as excinfo:
        ss.get_frame()
    assert (
        "Frames for orbits around custom bodies are not yet supported"
        in excinfo.exconly()
    )


@pytest.mark.parametrize(
    "body", [Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune]
)
def test_orbit_from_ephem_is_in_icrs_frame(body):
    ss = Orbit.from_body_ephem(body)

    assert ss.get_frame().is_equivalent_frame(ICRS())


def test_orbit_accepts_ecliptic_plane():
    r = [1e09, -4e09, -1e09] * u.km
    v = [5e00, -1e01, -4e00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v, plane=Planes.EARTH_ECLIPTIC)

    assert ss.get_frame().is_equivalent_frame(HeliocentricEclipticJ2000(obstime=J2000))


def test_orbit_represent_as_produces_correct_data():
    r = [1e09, -4e09, -1e09] * u.km
    v = [5e00, -1e01, -4e00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v)

    expected_result = CartesianRepresentation(
        *r, differentials=CartesianDifferential(*v)
    )

    result = ss.represent_as(CartesianRepresentation, CartesianDifferential)

    # We can't directly compare the objects, see
    # https://github.com/astropy/astropy/issues/7793
    assert (result.xyz == expected_result.xyz).all()
    assert (
        result.differentials["s"].d_xyz == expected_result.differentials["s"].d_xyz
    ).all()


def test_orbit_propagate_retains_plane():
    r = [1e09, -4e09, -1e09] * u.km
    v = [5e00, -1e01, -4e00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v, plane=Planes.EARTH_ECLIPTIC)

    orig_frame = ss.get_frame()

    final_ss = ss.propagate(1 * u.h)
    expected_frame = orig_frame.replicate_without_data(obstime=final_ss.epoch)

    assert final_ss.get_frame().is_equivalent_frame(expected_frame)


@pytest.mark.remote_data
def test_from_horizons_raise_valueerror():
    with pytest.raises(ValueError) as exep:
        Orbit.from_horizons(name="Dummy", attractor=Sun)
    assert (
        "ValueError: Unknown target (Dummy). Maybe try different id_type?"
        in exep.exconly()
    )


@pytest.mark.remote_data
def test_orbit_from_horizons_has_expected_elements():
    epoch = Time("2018-07-23", scale="tdb")
    # Orbit Parameters of Ceres
    # Taken from https://ssd.jpl.nasa.gov/horizons.cgi
    ss = Orbit.from_classical(
        Sun,
        2.76710759221651 * u.au,
        0.07554803091400027 * u.one,
        27.18502494739172 * u.deg,
        23.36913218336299 * u.deg,
        132.2919809219236 * u.deg,
        21.28957916690369 * u.deg,
        epoch,
    )
    ss1 = Orbit.from_horizons(name="Ceres", attractor=Sun, epoch=epoch)
    assert ss.pqw()[0].value.all() == ss1.pqw()[0].value.all()
    assert ss.r_a == ss1.r_a
    assert ss.a == ss1.a


@pytest.mark.remote_data
def test_plane_is_set_in_horizons():
    plane = Planes.EARTH_ECLIPTIC
    ss = Orbit.from_horizons(name="Ceres", attractor=Sun, plane=plane)
    assert ss.plane == plane


@pytest.mark.parametrize(
    "attractor,angular_velocity,expected_a,expected_period",
    [
        (
            Earth,
            (2 * np.pi / 23.9345) * u.rad / u.hour,
            42_164_205 * u.m,
            23.9345 * u.hour,
        ),
        (
            Mars,
            (2 * np.pi / 24.6228) * u.rad / u.hour,
            20_427_595 * u.m,
            24.6228 * u.hour,
        ),
    ],
)
def test_geostationary_creation_from_angular_velocity(
    attractor, angular_velocity, expected_a, expected_period
):
    ss = Orbit.geostationary(attractor=attractor, angular_velocity=angular_velocity)
    assert_quantity_allclose(ss.a, expected_a, rtol=1.0e-7)
    assert_quantity_allclose(ss.period, expected_period, rtol=1.0e-7)


@pytest.mark.parametrize(
    "attractor,period,expected_a",
    [
        (Earth, 23.9345 * u.hour, 42_164_205 * u.m),
        (Mars, 24.6228 * u.hour, 20_427_595 * u.m),
    ],
)
def test_geostationary_creation_from_period(attractor, period, expected_a):
    ss = Orbit.geostationary(attractor=attractor, period=period)
    assert_quantity_allclose(ss.a, expected_a, rtol=1.0e-7)
    assert_quantity_allclose(ss.period, period, rtol=1.0e-7)


@pytest.mark.parametrize(
    "attractor,period,hill_radius,expected_a",
    [
        (Earth, 23.9345 * u.hour, 0.01 * u.AU, 42_164_205 * u.m),
        (Mars, 24.6228 * u.hour, 1_000_000 * u.km, 20_427_595 * u.m),
    ],
)
def test_geostationary_creation_with_Hill_radius(
    attractor, period, hill_radius, expected_a
):
    ss = Orbit.geostationary(
        attractor=attractor, period=period, hill_radius=hill_radius
    )
    assert_quantity_allclose(ss.a, expected_a, rtol=1.0e-7)
    assert_quantity_allclose(ss.period, period, rtol=1.0e-7)


@pytest.mark.parametrize("attractor", [Earth, Mars])
def test_geostationary_input(attractor):
    with pytest.raises(ValueError) as excinfo:
        Orbit.geostationary(attractor=attractor)

    assert (
        "ValueError: At least one among angular_velocity or period must be passed"
        in excinfo.exconly()
    )


@pytest.mark.parametrize(
    "attractor,period,hill_radius", [(Venus, 243.025 * u.day, 1_000_000 * u.km)]
)
def test_geostationary_non_existence_condition(attractor, period, hill_radius):
    with pytest.raises(ValueError) as excinfo:
        Orbit.geostationary(attractor=attractor, period=period, hill_radius=hill_radius)

    assert (
        "Geostationary orbit for the given parameters doesn't exist"
        in excinfo.exconly()
    )


def test_heliosynchronous_orbit_enough_arguments():
    with pytest.raises(ValueError) as excinfo:
        Orbit.heliosynchronous(Earth, a=None, ecc=None, inc=None)

    assert (
        "At least two parameters of the set {a, ecc, inc} are required."
        in excinfo.exconly()
    )


def test_heliosynchronous_orbit_without_earth():
    with pytest.raises(NotImplementedError) as excinfo:
        Orbit.heliosynchronous(Mars, a=800 * u.km + Mars.R, ecc=0 * u.one)
    assert "Attractors other than Earth not supported yet" in excinfo.exconly()


def test_heliosynchronous_orbit_inc():
    # Vallado, example 11-2a
    expected_ecc = 0 * u.one
    expected_a = 800 * u.km + Earth.R
    expected_inc = 98.6 * u.deg
    ss0 = Orbit.heliosynchronous(Earth, a=expected_a, ecc=expected_ecc)

    assert_quantity_allclose(ss0.inc, expected_inc, rtol=1e-4)
    assert_quantity_allclose(ss0.a, expected_a)
    assert_quantity_allclose(ss0.ecc, expected_ecc)


def test_heliosynchronous_orbit_a():
    # Vallado, example 11-2b
    expected_ecc = 0.2 * u.one
    expected_inc = 98.6 * u.deg
    expected_a = 7346.846 * u.km
    ss0 = Orbit.heliosynchronous(Earth, ecc=expected_ecc, inc=expected_inc)

    assert_quantity_allclose(ss0.inc, expected_inc, rtol=1e-4)
    assert_quantity_allclose(ss0.a, expected_a, rtol=1e-5)
    assert_quantity_allclose(ss0.ecc, expected_ecc)


def test_perigee_and_apogee():
    expected_r_a = 500 * u.km
    expected_r_p = 300 * u.km
    a = (expected_r_a + expected_r_p) / 2
    ecc = expected_r_a / a - 1
    _a = 1.0 * u.deg  # Unused angle
    ss = Orbit.from_classical(Earth, a, ecc, _a, _a, _a, _a)
    assert_allclose(ss.r_a.to(u.km).value, expected_r_a.to(u.km).value)
    assert_allclose(ss.r_p.to(u.km).value, expected_r_p.to(u.km).value)


def test_expected_mean_anomaly():
    # Example from Curtis
    expected_mean_anomaly = 77.93 * u.deg

    attractor = Earth

    _a = 1.0 * u.deg  # Unused angle
    a = 15_300 * u.km
    ecc = 0.37255 * u.one
    nu = 120 * u.deg

    orbit = Orbit.from_classical(attractor, a, ecc, _a, _a, _a, nu)
    orbit_M = orbit.M

    assert_quantity_allclose(orbit_M.value, expected_mean_anomaly.value, rtol=1e-2)


def test_expected_angular_momentum():
    # Example from Curtis
    expected_ang_mag = 72472 * u.km ** 2

    attractor = Earth

    _a = 1.0 * u.deg  # Unused angle
    a = 15_300 * u.km
    ecc = 0.37255 * u.one
    nu = 120 * u.deg

    orbit = Orbit.from_classical(attractor, a, ecc, _a, _a, _a, nu)
    orbit_h_mag = orbit.h_mag

    assert_quantity_allclose(orbit_h_mag.value, expected_ang_mag.value, rtol=1e-2)


def test_expected_last_perifocal_passage():
    # Example from Curtis
    expected_t_p = 4077 * u.s

    attractor = Earth

    _a = 1.0 * u.deg  # Unused angle
    a = 15_300 * u.km
    ecc = 0.37255 * u.one
    nu = 120 * u.deg

    orbit = Orbit.from_classical(attractor, a, ecc, _a, _a, _a, nu)
    orbit_t_p = orbit.t_p

    assert_quantity_allclose(orbit_t_p.value, expected_t_p.value, rtol=1e-2)


def test_convert_from_rv_to_coe():
    # Data from Vallado, example 2.6
    attractor = Earth
    p = 11_067.790 * u.km
    ecc = 0.83285 * u.one
    inc = 87.87 * u.deg
    raan = 227.89 * u.deg
    argp = 53.38 * u.deg
    nu = 92.335 * u.deg
    expected_r = [6_525.344, 6_861.535, 6_449.125] * u.km
    expected_v = [4.902276, 5.533124, -1.975709] * u.km / u.s

    r, v = Orbit.from_classical(
        attractor, p / (1 - ecc ** 2), ecc, inc, raan, argp, nu
    ).rv()

    assert_quantity_allclose(r, expected_r, rtol=1e-5)
    assert_quantity_allclose(v, expected_v, rtol=1e-5)


def test_convert_from_coe_to_rv():
    # Data from Vallado, example 2.5
    attractor = Earth
    r = [6_524.384, 6_862.875, 6_448.296] * u.km
    v = [4.901327, 5.533756, -1.976341] * u.km / u.s

    expected_p = 11_067.79 * u.km
    expected_ecc = 0.832853 * u.one
    expected_inc = 87.870 * u.deg
    expected_raan = 227.89 * u.deg
    expected_argp = 53.38 * u.deg
    expected_nu = 92.335 * u.deg

    ss = Orbit.from_vectors(attractor, r, v)

    _, ecc, inc, raan, argp, nu = ss.classical()
    p = ss.p

    assert_quantity_allclose(p, expected_p, rtol=1e-4)
    assert_quantity_allclose(ecc, expected_ecc, rtol=1e-4)
    assert_quantity_allclose(inc, expected_inc, rtol=1e-4)
    assert_quantity_allclose(raan, expected_raan, rtol=1e-4)
    assert_quantity_allclose(argp, expected_argp, rtol=1e-4)
    assert_quantity_allclose(nu, expected_nu, rtol=1e-4)


def test_perifocal_points_to_perigee():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    ss = Orbit.from_classical(Sun, _d, _, _a, _a, _a, _a)
    p, _, _ = ss.pqw()
    assert_allclose(p, ss.e_vec / ss.ecc)


def test_arglat_within_range():
    r = [3_539.08827417, 5_310.19903462, 3_066.31301457] * u.km
    v = [-6.49780849, 3.24910291, 1.87521413] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)
    assert 0 * u.deg <= ss.arglat <= 360 * u.deg


def test_pqw_returns_dimensionless():
    r_0 = ([1, 0, 0] * u.au).to(u.km)  # type: ignore
    v_0 = ([0, 6, 0] * u.au / u.year).to(u.km / u.day)
    ss = Orbit.from_vectors(Sun, r_0, v_0)

    p, q, w = ss.pqw()

    assert p.unit == u.one
    assert q.unit == u.one
    assert w.unit == u.one


def test_from_coord_fails_if_no_time_differential():
    pos = [30_000, 0, 0] * u.km
    cartrep = CartesianRepresentation(*pos)

    # Method fails if coordinate instance doesn't contain a differential with
    # respect to time
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_coords(Earth, SkyCoord(cartrep))
    assert (
        "ValueError: Coordinate instance doesn't have a differential with respect to time"
        in excinfo.exconly()
    )


@pytest.mark.parametrize(
    "attractor", [Earth, Jupiter, Mars, Mercury, Neptune, Saturn, Sun, Uranus, Venus]
)
def test_orbit_creation_using_skycoord(attractor):
    vel = [0, 2, 0] * u.km / u.s
    cartdiff = CartesianDifferential(*vel)

    pos = [30_000, 0, 0] * u.km
    cartrep = CartesianRepresentation(*pos, differentials=cartdiff)

    coord = SkyCoord(cartrep, frame="icrs")
    o = Orbit.from_coords(attractor, coord)

    inertial_frame_at_body_centre = get_frame(
        attractor, Planes.EARTH_EQUATOR, obstime=coord.obstime
    )

    coord_transformed_to_irf = coord.transform_to(inertial_frame_at_body_centre)
    pos_transformed_to_irf = coord_transformed_to_irf.cartesian.xyz
    vel_transformed_to_irf = coord_transformed_to_irf.cartesian.differentials["s"].d_xyz

    assert (o.r == pos_transformed_to_irf).all()
    assert (o.v == vel_transformed_to_irf).all()


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "attractor", [Earth, Jupiter, Mars, Mercury, Neptune, Saturn, Sun, Uranus, Venus]
)
@pytest.mark.parametrize("frame", [ITRS, GCRS])
@pytest.mark.parametrize("obstime", [J2000, J2000_TDB])
def test_orbit_creation_using_frame_obj(attractor, frame, obstime):
    vel = [0, 2, 0] * u.km / u.s
    cartdiff = CartesianDifferential(*vel)

    pos = [30_000, 0, 0] * u.km
    cartrep = CartesianRepresentation(*pos, differentials=cartdiff)

    coord = frame(cartrep, obstime=obstime)
    o = Orbit.from_coords(attractor, coord)

    inertial_frame_at_body_centre = get_frame(
        attractor, Planes.EARTH_EQUATOR, obstime=coord.obstime
    )

    coord_transformed_to_irf = coord.transform_to(inertial_frame_at_body_centre)

    pos_transformed_to_irf = coord_transformed_to_irf.cartesian.xyz
    vel_transformed_to_irf = coord_transformed_to_irf.cartesian.differentials["s"].d_xyz

    assert_quantity_allclose(o.r, pos_transformed_to_irf, atol=1e-5 * u.km)
    assert_quantity_allclose(o.v, vel_transformed_to_irf, atol=1e-5 * u.km / u.s)


@pytest.mark.parametrize("obstime", [J2000, J2000_TDB])
def test_from_coord_fails_for_multiple_positions(obstime):
    cartdiff = CartesianDifferential(
        [[0, 1, 0], [-0.1, 0.9, 0]] * u.km / u.s, xyz_axis=1
    )
    cartrep = CartesianRepresentation(
        [[1, 0, 0], [0.9, 0.1, 0]] * u.km, differentials=cartdiff, xyz_axis=1
    )
    coords = GCRS(cartrep, representation_type=CartesianRepresentation, obstime=obstime)

    with pytest.raises(ValueError) as excinfo:
        Orbit.from_coords(Earth, coords)
    assert (
        "ValueError: Coordinate instance must represents exactly 1 position, found: 2"
        in excinfo.exconly()
    )


def test_from_coord_if_coord_is_not_of_shape_zero():
    pos = [0, 1, 0]
    vel = [1, 0, 0]
    cartdiff = CartesianDifferential([vel] * u.km / u.s, xyz_axis=1)
    cartrep = CartesianRepresentation([pos] * u.km, differentials=cartdiff, xyz_axis=1)
    coords = GCRS(cartrep, representation_type=CartesianRepresentation, obstime=J2000)

    ss = Orbit.from_coords(Earth, coords)

    assert_quantity_allclose(ss.r, pos * u.km, rtol=1e-5)
    assert_quantity_allclose(ss.v, vel * u.km / u.s, rtol=1e-5)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    "target_name", ["Ceres", "Vesta", "Eros"]
)  # Objects in both JPL SBDB and JPL Horizons
def test_from_sbdb_and_from_horizons_give_similar_results(target_name):
    ss_target = Orbit.from_sbdb(target_name)
    ss_classical = ss_target.classical()

    ss_ref = Orbit.from_horizons(
        name=target_name, attractor=Sun, plane=Planes.EARTH_ECLIPTIC
    )

    ss_ref = ss_ref.propagate_to_anomaly(
        ss_classical[5]
    )  # Catch reference orbit to same epoch
    ss_ref_class = ss_ref.classical()

    for test_elm, ref_elm in zip(ss_classical, ss_ref_class):
        assert_quantity_allclose(
            test_elm, ref_elm, rtol=1e-2
        )  # Maximum error of 1% (chosen arbitrarily)


@pytest.mark.remote_data
def test_from_sbdb_raise_valueerror():
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_sbdb(name="Halley")

    assert (
        str(excinfo.value)
        == "2 different objects found: \n2688 Halley (1982 HG1)\n1P/Halley\n"
    )


def test_orbit_wrong_dimensions_fails():
    bad_r = [[1000, 0, 0]] * u.km
    bad_v = [[[0, 10, 0]]] * u.km / u.s

    with pytest.raises(ValueError) as excinfo:
        Orbit.from_vectors(Earth, bad_r, bad_v)
    assert "ValueError: Vectors must have dimension 1, got 2 and 3" in excinfo.exconly()


@pytest.mark.remote_data
def test_orbit_change_attractor():
    Io = 501  # Id for Io moon
    ss_io = Orbit.from_horizons(Io, Sun, epoch=J2000, id_type="majorbody")
    ss_io = ss_io.change_attractor(Jupiter)
    assert Jupiter == ss_io.attractor


def test_orbit_change_attractor_returns_self():
    assert iss.change_attractor(iss.attractor) is iss


def test_orbit_change_attractor_out_of_SOI():
    ss = Orbit.from_vectors(
        Sun,
        r=[5.98967334e08, 4.09500684e08, 1.60955500e08] * u.km,
        v=[-13.30694373, 25.15256978, 11.59846936] * u.km / u.s,
        epoch=J2000,
    )

    with pytest.raises(ValueError) as excinfo:
        ss.change_attractor(Earth)
    assert "ValueError: Orbit is out of new attractor's SOI" in excinfo.exconly()


def test_orbit_change_attractor_force():
    ss = Orbit.from_vectors(
        Sun,
        r=[5.98967334e08, 4.09500684e08, 1.60955500e08] * u.km,
        v=[-13.30694373, 25.15256978, 11.59846936] * u.km / u.s,
        epoch=J2000,
    )

    ss_new_attractor = ss.change_attractor(Earth, force=True)
    assert ss_new_attractor.attractor == Earth


def test_orbit_change_attractor_unrelated_body():
    with pytest.raises(ValueError) as excinfo:
        iss.change_attractor(Mars)
    assert "ValueError: Cannot change to unrelated attractor" in excinfo.exconly()


def test_orbit_change_attractor_closed():
    with pytest.raises(ValueError) as excinfo:
        iss.change_attractor(Sun)
    assert (
        "ValueError: Orbit will never leave the SOI of its current attractor"
        in excinfo.exconly()
    )


def test_orbit_change_attractor_open():
    r = [-6_045, -3_490, 2_500] * u.km
    v = [-15.457, 6.618, 2.533] * u.km / u.s
    ss = Orbit.from_vectors(Earth, r, v)

    with pytest.warns(PatchedConicsWarning) as record:
        ss.change_attractor(Sun)

    w = record.pop(PatchedConicsWarning)
    assert "Leaving the SOI of the current attractor" in w.message.args[0]


@pytest.mark.parametrize(
    "expected_plane", [Planes.EARTH_ECLIPTIC, Planes.EARTH_EQUATOR]
)
def test_change_plane_sets_correct_plane(expected_plane):
    new_ss = iss.change_plane(expected_plane)

    assert new_ss.plane is expected_plane


def test_change_plane_twice_restores_original_data():
    new_ss = iss.change_plane(Planes.EARTH_ECLIPTIC).change_plane(iss.plane)

    assert_quantity_allclose(new_ss.r, iss.r)
    assert_quantity_allclose(new_ss.v, iss.v)


def test_time_to_anomaly():
    expected_tof = iss.period / 2
    iss_180 = iss.propagate_to_anomaly(180 * u.deg)
    tof = iss_180.time_to_anomaly(0 * u.deg)

    assert_quantity_allclose(tof, expected_tof)
