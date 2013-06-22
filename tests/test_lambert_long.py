# coding: utf-8
from poliastro.astiod import lambertuniv

r0 = [-1.09038587e8,  -1.04221731e8,   0.00000000e0]
rf = [-1.02959046e8,   3.05523274e7,   6.35994599e6]
tof = 5184000.0
k_Sun = 132712440018.0

va, vb, error = lambertuniv(r0, rf, 'L', 'N', tof, k_Sun)
