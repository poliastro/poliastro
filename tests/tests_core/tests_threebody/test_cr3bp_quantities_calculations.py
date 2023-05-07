from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
import pytest

from poliastro.bodies import Earth, Moon
from poliastro.core.threebody.cr3bp_quantities_calculations import (
    calculate_mu,
    calculate_tstar,
)


@pytest.mark.parametrize(
    "mu1, mu2, expected_mu",
    [
        (Earth.k, Moon.k, 1.215058560962404e-02 * u.one),
        # Earth-Moon expected mu (double click on systems): https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=halo&libr=1&branch=S
    ],
)
def test_calculate_mu(mu1, mu2, expected_mu):
    # Remember diffferent sources use different accuracy for GM values

    cr3bp_mu = calculate_mu(mu1, mu2)
    assert_quantity_allclose(cr3bp_mu, expected_mu, 1e-6)


@pytest.mark.parametrize(
    "mu1, mu2, lstar, expected_tstar",
    [
        (Earth.k, Moon.k, 389703.2648292776 * u.km, 382981.2891290545 * u.s),
        # Earth-Moon lstar (to compare tstar; double click on systems): https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=halo&libr=1&branch=S
        # This lstar is used instead of Moon.mean_a as https://ssd.jpl.nasa.gov/tools/periodic_orbits.html#/intro seems to use constants of slightly different accuracy
        # Earth-Moon expected tstar (double click on systems): https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=halo&libr=1&branch=S
    ],
)
def test_calculate_tstar(mu1, mu2, lstar, expected_tstar):
    # Earth-Moon mu: https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=halo&libr=1&branch=S
    # Remember diffferent sources use different accuracy for GM values

    cr3bp_tstar = calculate_tstar(mu1, mu2, lstar)
    assert_quantity_allclose(cr3bp_tstar, expected_tstar, 1e-6)
