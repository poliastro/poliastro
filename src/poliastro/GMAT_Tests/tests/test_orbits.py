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
    report = pd.read_csv(os.path.join(report_path, "Cartesian_orbit.csv"))

    initial_vector = report.iloc[0].values
    final_vector = report.iloc[-1].values
    r_initial = initial_vector[1:4] * u.km
    v_initial = initial_vector[4:] * u.km / u.s
    # tof is 12000s
    tof = final_vector[0] * u.s
    ss = Orbit.from_vectors(Earth, r_initial, v_initial)
    position_final = propagate(ss, time.TimeDelta(tof))
    a = np.array(
        [
            position_final.x.value,
            position_final.y.value,
            position_final.z.value,
            position_final.v_x.value,
            position_final.v_y.value,
            position_final.v_z.value,
        ]
    ).flatten()
    for i, j in zip(a, final_vector[1:]):
        assert_quantity_allclose(i, j)
