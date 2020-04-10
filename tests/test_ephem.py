from unittest import mock

import pytest
from astropy import units as u
from astropy.coordinates import (
    ICRS,
    BarycentricMeanEcliptic,
    CartesianDifferential,
    CartesianRepresentation,
)
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from poliastro.bodies import Earth, Venus
from poliastro.ephem import Ephem, InterpolationMethods
from poliastro.frames import Planes

AVAILABLE_INTERPOLATION_METHODS = InterpolationMethods.__members__.values()
AVAILABLE_PLANES = Planes.__members__.values()


def assert_coordinates_allclose(actual, desired, rtol=1e-7, atol_scale=None, **kwargs):
    if atol_scale is None:
        atol_scale = 0

    assert_quantity_allclose(
        actual.xyz, desired.xyz, rtol, atol=atol_scale * desired.xyz.unit, **kwargs
    )
    if "s" in desired.differentials:
        assert_quantity_allclose(
            actual.differentials["s"].d_xyz,
            desired.differentials["s"].d_xyz,
            rtol=rtol,
            atol=atol_scale * desired.differentials["s"].d_xyz.unit,
            **kwargs,
        )


@pytest.fixture
def epochs():
    return Time(
        [
            "2020-03-01 12:00:00",
            "2020-03-02 12:00:00",
            "2020-03-03 12:00:00",
            "2020-03-04 12:00:00",
        ],
        scale="tdb",
    )


@pytest.fixture
def coordinates():
    return CartesianRepresentation(
        [(1, 0, 0), (0.9, 0.1, 0), (0.8, 0.2, 0), (0.7, 0.3, 0)] * u.au,
        xyz_axis=1,
        differentials=CartesianDifferential(
            [(0, 1, 0), (-0.1, 0.9, 0), (-0.2, 0.8, 0), (-0.3, 0.7, 0)]
            * (u.au / u.day),
            xyz_axis=1,
        ),
    )


@pytest.mark.parametrize("plane", AVAILABLE_PLANES)
def test_ephem_has_given_plane(plane):
    ephem = Ephem(None, None, plane)

    assert ephem.plane is plane


@pytest.mark.parametrize("method", AVAILABLE_INTERPOLATION_METHODS)
def test_ephem_sample_no_arguments_returns_exactly_same_input(
    epochs, coordinates, method
):
    ephem = Ephem(coordinates, epochs, None)

    result_coordinates = ephem.sample(method=method)

    # Exactly the same
    assert result_coordinates == coordinates


@pytest.mark.parametrize("method", AVAILABLE_INTERPOLATION_METHODS)
def test_ephem_sample_same_epochs_returns_same_input(epochs, coordinates, method):
    ephem = Ephem(coordinates, epochs, None)

    result_coordinates = ephem.sample(epochs, method=method)

    # TODO: Should it return exactly the same?
    assert_coordinates_allclose(result_coordinates, coordinates, atol_scale=1e-17)


@pytest.mark.parametrize("method", AVAILABLE_INTERPOLATION_METHODS)
def test_ephem_sample_existing_epochs_returns_corresponding_input(
    epochs, coordinates, method
):
    ephem = Ephem(coordinates, epochs, None)

    result_coordinates = ephem.sample(epochs[::2], method=method)

    # Exactly the same
    assert_coordinates_allclose(result_coordinates, coordinates[::2], atol_scale=1e-17)


def test_rv_no_parameters_returns_input_vectors(coordinates, epochs):
    unused_plane = Planes.EARTH_EQUATOR
    ephem = Ephem(coordinates, epochs, unused_plane)

    expected_r = coordinates.get_xyz(xyz_axis=1)
    expected_v = coordinates.differentials["s"].get_d_xyz(xyz_axis=1)

    r, v = ephem.rv()

    assert_quantity_allclose(r, expected_r)
    assert_quantity_allclose(v, expected_v)


def test_rv_scalar_epoch_returns_scalar_vectors(coordinates, epochs):
    unused_plane = Planes.EARTH_EQUATOR
    ephem = Ephem(coordinates, epochs, unused_plane)

    expected_r = coordinates.get_xyz(xyz_axis=1)[0]
    expected_v = coordinates.differentials["s"].get_d_xyz(xyz_axis=1)[0]

    r, v = ephem.rv(epochs[0])

    assert_quantity_allclose(r, expected_r)
    assert_quantity_allclose(v, expected_v)


@pytest.mark.parametrize("method", AVAILABLE_INTERPOLATION_METHODS)
@pytest.mark.parametrize(
    "plane, FrameClass, rtol",
    [
        (Planes.EARTH_EQUATOR, ICRS, 1e-7),
        (Planes.EARTH_ECLIPTIC, BarycentricMeanEcliptic, 1e-5),
    ],
)
def test_ephem_from_body_has_expected_properties(method, plane, FrameClass, rtol):
    epochs = Time(
        ["2020-03-01 12:00:00", "2020-03-17 00:00:00.000", "2020-04-01 12:00:00.000"],
        scale="tdb",
    )
    equatorial_coordinates = CartesianRepresentation(
        [
            (-1.40892271e08, 45067626.83900666, 19543510.68386639),
            (-1.4925067e08, 9130104.71634121, 3964948.59999307),
            (-1.46952333e08, -27413113.24215863, -11875983.21773582),
        ]
        * u.km,
        xyz_axis=1,
        differentials=CartesianDifferential(
            [
                (-10.14262131, -25.96929533, -11.25810932),
                (-2.28639444, -27.3906416, -11.87218591),
                (5.67814544, -26.84316701, -11.63720607),
            ]
            * (u.km / u.s),
            xyz_axis=1,
        ),
    )

    expected_coordinates = (
        ICRS(equatorial_coordinates)
        .transform_to(FrameClass)
        .represent_as(CartesianRepresentation, CartesianDifferential)
    )

    earth = Ephem.from_body(Earth, epochs, plane)
    coordinates = earth.sample(method=method)

    assert earth.epochs is epochs
    assert_coordinates_allclose(coordinates, expected_coordinates, rtol=rtol)


def test_from_body_scalar_epoch_uses_reshaped_epochs():
    expected_epochs = Time(["2020-03-01 12:00:00"], scale="tdb")
    epochs = expected_epochs[0]

    unused_plane = Planes.EARTH_EQUATOR
    ephem = Ephem.from_body(Earth, epochs, plane=unused_plane)

    assert ephem.epochs == expected_epochs


@mock.patch("poliastro.ephem.Horizons")
@pytest.mark.parametrize("body,location_str", [(Earth, "500@399"), (Venus, "500@299")])
@pytest.mark.parametrize(
    "plane,refplane_str",
    [(Planes.EARTH_EQUATOR, "earth"), (Planes.EARTH_ECLIPTIC, "ecliptic")],
)
def test_ephem_from_horizons_calls_horizons_with_correct_parameters(
    horizons_mock, body, location_str, plane, refplane_str
):
    unused_name = "Strange Object"
    unused_id_type = "id_type"
    epochs = Time(["2020-03-01 12:00:00"], scale="tdb")

    horizons_mock().vectors.return_value = {
        "x": [1] * u.au,
        "y": [0] * u.au,
        "z": [0] * u.au,
        "vx": [0] * (u.au / u.day),
        "vy": [1] * (u.au / u.day),
        "vz": [0] * (u.au / u.day),
    }
    expected_coordinates = CartesianRepresentation(
        [(1, 0, 0)] * u.au,
        xyz_axis=1,
        differentials=CartesianDifferential([(0, 1, 0)] * (u.au / u.day), xyz_axis=1),
    )

    ephem = Ephem.from_horizons(unused_name, body, epochs, plane, unused_id_type)

    horizons_mock.assert_called_with(
        id=unused_name, location=location_str, epochs=epochs.jd, id_type=unused_id_type
    )
    horizons_mock().vectors.assert_called_once_with(refplane=refplane_str)

    coordinates = ephem.sample()

    assert_coordinates_allclose(coordinates, expected_coordinates)


@mock.patch("poliastro.ephem.Horizons")
def test_from_horizons_scalar_epoch_uses_reshaped_epochs(horizons_mock):
    unused_name = "Strange Object"
    unused_id_type = "id_type"
    unused_plane = Planes.EARTH_EQUATOR
    unused_location_str = "500@399"
    unused_attractor = Earth

    expected_epochs = Time(["2020-03-01 12:00:00"], scale="tdb")
    epochs = expected_epochs[0]

    horizons_mock().vectors.return_value = {
        "x": [1] * u.au,
        "y": [0] * u.au,
        "z": [0] * u.au,
        "vx": [0] * (u.au / u.day),
        "vy": [1] * (u.au / u.day),
        "vz": [0] * (u.au / u.day),
    }

    Ephem.from_horizons(
        unused_name,
        epochs,
        attractor=unused_attractor,
        plane=unused_plane,
        id_type=unused_id_type,
    )

    horizons_mock.assert_called_with(
        id=unused_name,
        location=unused_location_str,
        epochs=expected_epochs.jd,
        id_type=unused_id_type,
    )
