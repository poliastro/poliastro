import astropy.units as u
from astropy import time
from astropy.tests.helper import assert_quantity_allclose

from poliastro import coordinates, bodies
from poliastro.constants import J2000


# Note that function are tested using astropy current builtin ephemeris.
# Horizons uses JPL ephemeris DE431, so expected values are hardcoded,
# instead of being obtained using Horizons.
def test_body_centered_to_icrs_transformation():

    vexpress_r_venus = [-2.707041558060933E+03, 1.112962479175306E+04, -3.792944408664889E+04] * u.km
    vexpress_v_venus = [-2.045118200275925E-01, 7.978578482960554E-01, 2.664944903217139E+00] * u.km / u.s

    expected_r = [-3.47202219448080286E+07, 9.16853879708216339E+07, 4.34117810525591150E+07] * u.km
    expected_v = [-3.34053728321152121E+01, -1.16604776013667291E+01, -2.39943678872506838E-01] * u.km / u.s

    r, v = coordinates.body_centered_to_icrs(vexpress_r_venus, vexpress_v_venus, bodies.Venus,
                                             time.Time("2014-08-23 00:00", scale='tdb'))

    assert_quantity_allclose(r, expected_r)
    assert_quantity_allclose(v, expected_v)


def test_icrs_to_body_centered_transformation():
    vexpress_r_icrs = [-3.472125578094885E+07, 9.168528034176737E+07, 4.341160627674723E+07] * u.km
    vexpress_v_icrs = [-3.340574196483147E+01, -1.165974037637970E+01, -2.395829145441408E-01] * u.km / u.s

    expected_r = [-3.74486105008138566e+03, 1.10085874027602295e+04, -3.80681106516677464e+04] * u.km
    expected_v = [-2.04845025352488774e-01, 7.98692896032012989e-01, 2.66498465286454023e+00] * u.km / u.s

    r, v = coordinates.icrs_to_body_centered(vexpress_r_icrs, vexpress_v_icrs, bodies.Venus,
                                             time.Time("2014-08-23 00:00", scale='tdb'))

    assert_quantity_allclose(r, expected_r)
    assert_quantity_allclose(v, expected_v)
