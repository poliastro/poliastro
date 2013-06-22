# Test Lambert wrapped subroutines
# TODO: This interface without k is very, very stupid and leads to errors:
# change

import numpy as np
import poliastro.astiod

r0 = np.array([15945.15, 0, 0])
r = np.array([12214.83899, 10249.46731, 0])
dt = 76.0 * 60

v0, v, error = poliastro.astiod.lambertuniv(r0, r, 's', 'n', dt)
print(v0, v)

v0, v, error = poliastro.astiod.lambertbattin(r0, r, 's', 'n', dt)
print(v0, v)
