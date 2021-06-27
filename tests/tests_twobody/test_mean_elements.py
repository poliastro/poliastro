import pytest

from poliastro.twobody.mean_elements import get_mean_elements


def test_get_mean_elements_raises_error_if_invalid_body():
    body = "Sun"

    with pytest.raises(ValueError) as excinfo:
        get_mean_elements(body)
    assert f"The input body '{body}' is invalid." in excinfo.exconly()
