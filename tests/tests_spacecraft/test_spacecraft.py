from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
import numpy as np

from poliastro.spacecraft import Spacecraft


def test_spacecraft_init():
    C_D = 2.2 * u.one  # Dimensionless (any value would do)
    A = ((np.pi / 4.0) * (u.m**2)).to(u.km**2)
    m = 100 * u.kg
    spacecraft = Spacecraft(A, C_D, m)
    assert isinstance(spacecraft, Spacecraft)


def test_balistic_coefficient():
    C_D = 2.2 * u.one  # Dimensionless (any value would do)
    A = ((np.pi / 4.0) * (u.m**2)).to(u.km**2)
    m = 100 * u.kg
    spacecraft = Spacecraft(A, C_D, m)
    assert_quantity_allclose(
        spacecraft.ballistic_coefficient.to_value(),
        1.7278759594743866e-08,
        rtol=1e-10,
    )
