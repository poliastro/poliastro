from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
import pytest

from poliastro.core.propagation import (
    danby_coe,
    gooding_coe,
    markley_coe,
    mikkola_coe,
    pimienta_coe,
)
from poliastro.core.propagation.farnocchia import farnocchia_coe
from poliastro.examples import iss


@pytest.mark.parametrize(
    "propagator_coe",
    [
        danby_coe,
        markley_coe,
        pimienta_coe,
        mikkola_coe,
        farnocchia_coe,
        gooding_coe,
    ],
)
def test_propagate_with_coe(propagator_coe):

    period = iss.period
    a, ecc, inc, raan, argp, nu = iss.classical()
    p = a * (1 - ecc**2)

    # Delete the units
    p = p.to_value(u.km)
    ecc = ecc.value
    period = period.to_value(u.s)
    inc = inc.to_value(u.rad)
    raan = raan.to_value(u.rad)
    argp = argp.to_value(u.rad)
    nu = nu.to_value(u.rad)
    k = iss.attractor.k.to_value(u.km**3 / u.s**2)

    nu_final = propagator_coe(k, p, ecc, inc, raan, argp, nu, period)

    assert_quantity_allclose(nu_final, nu)
