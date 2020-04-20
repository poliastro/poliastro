from astropy.tests.helper import assert_quantity_allclose
from numpy.linalg import norm

from poliastro import bodies, coordinates
from poliastro.examples import molniya
from poliastro.twobody.orbit import Orbit


def test_inertial_body_centered_to_pqw():
    molniya_r_peri, molniya_v_peri = coordinates.inertial_body_centered_to_pqw(
        molniya.r, molniya.v, bodies.Earth
    )

    molniya_peri = Orbit.from_vectors(
        bodies.Earth, molniya_r_peri, molniya_v_peri, molniya.epoch
    )

    assert_quantity_allclose(molniya_peri.e_vec[-2:], [0, 0], atol=1e-12)
    assert_quantity_allclose(norm(molniya_peri.e_vec), norm(molniya.e_vec))
