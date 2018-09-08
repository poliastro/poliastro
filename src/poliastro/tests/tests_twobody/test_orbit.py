import pytest

import pickle

from numpy.testing import assert_allclose, assert_array_equal

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import CartesianRepresentation, CartesianDifferential

from astropy.time import Time

from poliastro.bodies import (
    Body,
    Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto,
)
from poliastro.frames import (
    Planes,
    ICRS,
    HCRS, MercuryICRS, VenusICRS, GCRS, MarsICRS, JupiterICRS, SaturnICRS, UranusICRS, NeptuneICRS, PlutoICRS,
    HeliocentricEclipticJ2000
)
from poliastro.twobody import Orbit
from poliastro.twobody.orbit import TimeScaleWarning
from poliastro.constants import J2000


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
    assert ("UnitsError: Argument 'nu' to function 'from_classical' must be in units convertible to 'rad'."
            in excinfo.exconly())


def test_state_raises_unitserror_if_rv_units_are_wrong():
    _d = [1.0, 0.0, 0.0] * u.AU
    wrong_v = [0.0, 1.0e-6, 0.0] * u.AU
    with pytest.raises(u.UnitsError) as excinfo:
        Orbit.from_vectors(Sun, _d, wrong_v)
    assert ("UnitsError: Argument 'v' to function 'from_vectors' must be in units convertible to 'm / s'."
            in excinfo.exconly())


def test_parabolic_elements_fail_early():
    attractor = Earth
    ecc = 1.0 * u.one
    _d = 1.0 * u.AU  # Unused distance
    _a = 1.0 * u.deg  # Unused angle
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_classical(attractor, _d, ecc, _a, _a, _a, _a)
    assert ("ValueError: For parabolic orbits use Orbit.parabolic instead" in excinfo.exconly())


def test_bad_inclination_raises_exception():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    bad_inc = 200 * u.deg
    _body = Sun  # Unused body
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_classical(_body, _d, _, bad_inc, _a, _a, _a)
    assert ("ValueError: Inclination must be between 0 and 180 degrees" in excinfo.exconly())


def test_bad_hyperbolic_raises_exception():
    bad_a = 1.0 * u.AU
    ecc = 1.5 * u.one
    _a = 1.0 * u.deg  # Unused angle
    _inc = 100 * u.deg  # Unused inclination
    _body = Sun  # Unused body
    with pytest.raises(ValueError) as excinfo:
        Orbit.from_classical(_body, bad_a, ecc, _inc, _a, _a, _a)
    assert ("Hyperbolic orbits have negative semimajor axis" in excinfo.exconly())


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


def test_orbit_representation():
    ss = Orbit.circular(Earth, 600 * u.km, 20 * u.deg)
    expected_str = "6978 x 6978 km x 20.0 deg (GCRS) orbit around Earth (\u2641)"

    assert str(ss) == repr(ss) == expected_str


def test_sample_numpoints():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    ss = Orbit.from_classical(_body, _d, _, _a, _a, _a, _a)
    positions = ss.sample(values=50)
    assert len(positions) == 50


def test_sample_with_time_value():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    ss = Orbit.from_classical(_body, _d, _, _a, _a, _a, _a)

    expected_r = [ss.r]
    positions = ss.sample(values=ss.nu + [360] * u.deg)
    r = positions.data.xyz.transpose()

    assert_quantity_allclose(r, expected_r, rtol=1.e-7)


def test_sample_with_nu_value():
    _d = 1.0 * u.AU  # Unused distance
    _ = 0.5 * u.one  # Unused dimensionless value
    _a = 1.0 * u.deg  # Unused angle
    _body = Sun  # Unused body
    ss = Orbit.from_classical(_body, _d, _, _a, _a, _a, _a)

    expected_r = [ss.r]
    positions = ss.sample(values=ss.nu + [360] * u.deg)
    r = positions.data.xyz.transpose()

    assert_quantity_allclose(r, expected_r, rtol=1.e-7)


def test_hyperbolic_nu_value_check():
    # A custom hyperbolic orbit
    r = [1.197659243752796E+09, -4.443716685978071E+09, -1.747610548576734E+09] * u.km
    v = [5.540549267188614E+00, -1.251544669134140E+01, -4.848892572767733E+00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v, Time('2015-07-14 07:59', scale='tdb'))

    positions = ss.sample(100)

    assert isinstance(positions, HCRS)
    assert len(positions) == 100


def test_hyperbolic_modulus_wrapped_nu():
    ss = Orbit.from_vectors(
        Sun,
        [-9.77441841e+07, 1.01000539e+08, 4.37584668e+07] * u.km,
        [23.75936985, -43.09599568, -8.7084724] * u.km / u.s,
    )
    num_values = 3

    positions = ss.sample(num_values)

    assert_quantity_allclose(positions[0].data.xyz, ss.r)


def test_orbit_is_pickable():
    # A custom hyperbolic orbit
    r = [1.197659243752796E+09, -4.443716685978071E+09, -1.747610548576734E+09] * u.km
    v = [5.540549267188614E+00, -1.251544669134140E+01, -4.848892572767733E+00] * u.km / u.s
    epoch = Time('2015-07-14 07:59', scale='tdb')

    ss = Orbit.from_vectors(Sun, r, v, epoch)

    pickled = pickle.dumps(ss)
    ss_result = pickle.loads(pickled)

    assert_array_equal(ss.r, ss_result.r)
    assert_array_equal(ss.v, ss_result.v)
    assert ss_result.epoch == ss.epoch


@pytest.mark.parametrize("attractor, expected_frame_class", [
    (Sun, HCRS),
    (Mercury, MercuryICRS),
    (Venus, VenusICRS),
    pytest.param(Earth, GCRS, marks=pytest.mark.xfail),  # See https://github.com/astropy/astropy/issues/7793
    (Mars, MarsICRS),
    (Jupiter, JupiterICRS),
    (Saturn, SaturnICRS),
    (Uranus, UranusICRS),
    (Neptune, NeptuneICRS),
    (Pluto, PlutoICRS),
])
def test_orbit_has_proper_frame(attractor, expected_frame_class):
    # Dummy data
    r = [1E+09, -4E+09, -1E+09] * u.km
    v = [5E+00, -1E+01, -4E+00] * u.km / u.s
    epoch = Time('2015-07-14 07:59', scale='tdb')

    ss = Orbit.from_vectors(attractor, r, v, epoch)

    assert ss.frame.is_equivalent_frame(expected_frame_class(obstime=epoch))
    assert ss.frame.obstime == epoch


def test_orbit_from_custom_body_raises_error_when_asked_frame():
    attractor = Body(Sun, 1 * u.km ** 3 / u.s ** 2, "_DummyPlanet")

    r = [1E+09, -4E+09, -1E+09] * u.km
    v = [5E+00, -1E+01, -4E+00] * u.km / u.s

    ss = Orbit.from_vectors(attractor, r, v)

    with pytest.raises(NotImplementedError) as excinfo:
        ss.frame
    assert ("Frames for orbits around custom bodies are not yet supported"
            in excinfo.exconly())


@pytest.mark.parametrize("body", [
    Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
])
def test_orbit_from_ephem_is_in_icrs_frame(body):
    ss = Orbit.from_body_ephem(body)

    assert ss.frame.is_equivalent_frame(ICRS())


def test_orbit_accepts_ecliptic_plane():
    r = [1E+09, -4E+09, -1E+09] * u.km
    v = [5E+00, -1E+01, -4E+00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v, plane=Planes.EARTH_ECLIPTIC)

    assert ss.frame.is_equivalent_frame(HeliocentricEclipticJ2000(obstime=J2000))


def test_orbit_represent_as_produces_correct_data():
    r = [1E+09, -4E+09, -1E+09] * u.km
    v = [5E+00, -1E+01, -4E+00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v)

    expected_result = CartesianRepresentation(
        *r, differentials=CartesianDifferential(*v)
    )

    result = ss.represent_as(CartesianRepresentation)

    # We can't directly compare the objects, see https://github.com/astropy/astropy/issues/7793
    assert (result.xyz == expected_result.xyz).all()
    assert (result.differentials['s'].d_xyz == expected_result.differentials['s'].d_xyz).all()


@pytest.mark.parametrize("value", [1 * u.h, 10 * u.deg])
def test_orbit_propagate_retains_plane(value):
    r = [1E+09, -4E+09, -1E+09] * u.km
    v = [5E+00, -1E+01, -4E+00] * u.km / u.s

    ss = Orbit.from_vectors(Sun, r, v, plane=Planes.EARTH_ECLIPTIC)

    orig_frame = ss.frame

    final_ss = ss.propagate(1 * u.h)
    expected_frame = orig_frame.replicate_without_data(obstime=final_ss.epoch)

    assert final_ss.frame.is_equivalent_frame(expected_frame)
