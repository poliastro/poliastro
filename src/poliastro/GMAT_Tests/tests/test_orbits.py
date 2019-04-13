import os

import numpy as np
import pandas as pd
from astropy import time, units as u
from astropy.tests.helper import assert_quantity_allclose

from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import propagate


def test_from_vectors():
    report_path = os.path.join(os.getcwd(), "src/poliastro/GMAT_Tests/Reports")
    with open(os.path.join(report_path, "Cartesian_orbit.txt")) as f:
        report = f.readlines()
    # report of the form [elapsed seconds X, Y, Z, VX, VY, VZ]
    initial_vector_GMAT = [float(i) for i in report[1].strip("\n").split(" ")]
    final_vector_GMAT = [float(i) for i in report[-1].strip("\n").split(" ")]
    r_initial = initial_vector_GMAT[1:4] * u.km
    v_initial = initial_vector_GMAT[4:] * u.km / u.s
    # tof is 12000s
    tof = final_vector_GMAT[0] * u.s
    ss = Orbit.from_vectors(Earth, r_initial, v_initial)
    position_final = propagate(ss, time.TimeDelta(tof))
    final_vector_poliastro = np.array(
        [
            position_final.x.value,
            position_final.y.value,
            position_final.z.value,
            position_final.v_x.value,
            position_final.v_y.value,
            position_final.v_z.value,
        ]
    ).flatten()
    for i, j in zip(final_vector_poliastro, final_vector_GMAT[1:]):
        assert_quantity_allclose(i, j)


def test_from_classical():
    report_path = os.path.join(os.getcwd(), "src/poliastro/GMAT_Tests/Reports")
    with open(os.path.join(report_path, "Keplerian_orbit.txt")) as f:
        report = f.readlines()
    # report of the form [elapsed seconds X, Y, Z, VX, VY, VZ, SMA, ECC, INC, RAAN, ARPG, NU]
    initial_vector_GMAT = [float(i) for i in report[1].strip("\n").split(" ")]
    final_vector_GMAT = [float(i) for i in report[-1].strip("\n").split(" ")]
    r_initial = initial_vector_GMAT[1:4] * u.km
    v_initial = initial_vector_GMAT[4:7] * u.km / u.s
    # tof is 12000s
    print(r_initial, v_initial)
    tof = final_vector_GMAT[0] * u.s
    a, ecc, inc, raan, arpg, nu = initial_vector_GMAT[7:]
    ss = Orbit.from_classical(
        Earth,
        a * u.km,
        ecc * u.one,
        inc * u.deg,
        raan * u.deg,
        arpg * u.deg,
        nu * u.deg,
    )
    position_final = propagate(ss, time.TimeDelta(tof))
    final_vector_poliastro = np.array(
        [
            position_final.x.value,
            position_final.y.value,
            position_final.z.value,
            position_final.v_x.value,
            position_final.v_y.value,
            position_final.v_z.value,
        ]
    ).flatten()
    for i, j in zip(final_vector_poliastro, final_vector_GMAT[1:]):
        assert_quantity_allclose(i, j)
