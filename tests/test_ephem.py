from astropy import units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from poliastro.bodies import Earth
from poliastro.ephem import Ephem


def assert_coordinates_allclose(actual, desired, *args, **kwargs):
    assert_quantity_allclose(actual.xyz, desired.xyz, *args, **kwargs)
    if "s" in desired.differentials:
        assert_quantity_allclose(
            actual.differentials["s"].d_xyz,
            desired.differentials["s"].d_xyz,
            *args,
            **kwargs,
        )


def test_ephem_from_body_has_expected_properties():
    epochs = Time(
        ["2020-03-01 12:00:00", "2020-03-17 00:00:00.000", "2020-04-01 12:00:00.000"],
        scale="tdb",
    )
    expected_coordinates = CartesianRepresentation(
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

    earth = Ephem.from_body(Earth, epochs)
    coordinates = earth.sample()

    assert earth.epochs is epochs
    assert_coordinates_allclose(coordinates, expected_coordinates)
