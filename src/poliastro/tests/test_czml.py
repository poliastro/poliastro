import pytest
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from astropy.time import TimeDelta
from numpy import insert as np_insert
from numpy.testing import assert_allclose

from poliastro.czml.extract_czml import CZMLExtractor
from poliastro.examples import iss, molniya
from poliastro.twobody.propagation import propagate


def test_czml_add_orbit():
    start_epoch = iss.epoch
    end_epoch = iss.epoch + molniya.period

    sample_points = 10

    extractor = CZMLExtractor(start_epoch, end_epoch, sample_points)

    extractor.add_orbit(
        molniya, label_text="Molniya", label_fill_color=[125, 80, 120, 255]
    )
    extractor.add_orbit(iss, label_text="ISS", path_show=False)

    cords_iss = extractor.czml[1]["position"]["cartesian"]

    h_iss = (end_epoch - iss.epoch).to(u.second) / sample_points

    for i in range(sample_points):
        position_iss_test = propagate(iss, TimeDelta(i * h_iss))
        cords_iss_test = (
            position_iss_test.represent_as(CartesianRepresentation)
            .xyz.to(u.meter)
            .value
        )
        cords_iss_test = np_insert(cords_iss_test, 0, h_iss.value * i, axis=0)

        for j in range(4):
            assert_allclose(cords_iss[4 * i + j], cords_iss_test[j], rtol=1e-5)

    # Test label params and that the values between the two objects are not overwritten
    assert extractor.czml[0]["label"]["text"] == "Molniya"
    assert extractor.czml[1]["label"]["text"] == "ISS"
    assert extractor.czml[0]["label"]["fillColor"]["rgba"] == [125, 80, 120, 255]

    assert extractor.czml[1]["path"]["show"]["boolean"] is False


def test_czml_invalid_orbit_epoch_error():
    start_epoch = molniya.epoch
    end_epoch = molniya.epoch + molniya.period

    extractor = CZMLExtractor(start_epoch, end_epoch, 10)

    with pytest.raises(ValueError) as excinfo:
        extractor.add_orbit(iss, label_text="ISS", path_show=False)
    assert (
        "ValueError: The orbit's epoch cannot exceed the constructor's ending epoch"
        in excinfo.exconly()
    )
