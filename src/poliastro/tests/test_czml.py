import pytest
from numpy.testing import assert_allclose

from poliastro.czml.extract_czml import CZMLExtractor
from poliastro.examples import iss, molniya


def test_czml_add_orbit():
    start_epoch = iss.epoch
    end_epoch = iss.epoch + molniya.period

    Extractor = CZMLExtractor(start_epoch, end_epoch, 10)

    Extractor.add_orbit(
        molniya, label_text="Molniya", label_fill_color=[125, 80, 120, 255]
    )
    Extractor.add_orbit(iss, label_text="ISS", path_show=False)

    cords_molniya = Extractor.czml[0]["position"]["cartesian"][4 * 5 : 4 * 6]
    cords_molniya_test = [
        21587.554141072753,
        17029383.601313736,
        5939060.590166127,
        11860029.944171576,
    ]

    cords_iss = Extractor.czml[1]["position"]["cartesian"][4 * 5 : 4 * 6]
    cords_iss_test = [
        21587.554141072753,
        -3674993.3280015206,
        -4315653.729708829,
        3705413.324941563,
    ]

    # Test orbital coordinates
    for i in range(4):
        assert_allclose(cords_molniya[i], cords_molniya_test[i], rtol=1e-5)
        assert_allclose(cords_iss[i], cords_iss_test[i], rtol=1e-5)

    # Test label params and that the values between the two objects are not overwritten
    assert Extractor.czml[0]["label"]["text"] == "Molniya"
    assert Extractor.czml[1]["label"]["text"] == "ISS"
    assert Extractor.czml[0]["label"]["fillColor"]["rgba"] == [125, 80, 120, 255]

    assert Extractor.czml[1]["path"]["show"]["boolean"] is False


def test_czml_invalid_orbit_epoch_error():
    start_epoch = molniya.epoch
    end_epoch = molniya.epoch + molniya.period

    Extractor = CZMLExtractor(start_epoch, end_epoch, 10)

    with pytest.raises(ValueError) as excinfo:
        Extractor.add_orbit(iss, label_text="ISS", path_show=False)
    assert (
        "ValueError: The orbit's epoch cannot exceed the constructor's ending epoch"
        in excinfo.exconly()
    )
