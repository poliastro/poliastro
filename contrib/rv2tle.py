"""
Author: jtegedor

Python implementation of RV2TLE program for computation of mean orbital elements
from state vector, for an Earth-orbitting satellite

Original C++ implementation available in http://sat.belastro.net/satelliteorbitdetermination.com/.
GitHub repository with MIT license here: https://github.com/interplanetarychris/scottcampbell-satfit

Variables names have been kept similar to the original implementation, for easier comparison
"""

import math

import numpy as np
from astropy import coordinates as coord, units as u
from astropy.coordinates import TEME
from astropy.time import Time
from sgp4.api import SGP4_ERRORS, WGS84, Satrec

from poliastro.bodies import Earth
from poliastro.examples import iss
from poliastro.twobody import Orbit, angles


def unitv(v):
    """
    Compute unitary vector of numpy array v
    """
    return v / np.linalg.norm(v)


def acose(x):
    """
    Custom implementation of numpy.acose
    Returns either 0 or math.pi if input value out of [-1, 1]
    """
    result = 0.0
    if x > 1:
        result = 0
    elif x < -1:
        result = math.pi
    else:
        result = math.acos(x)
    return result


def fmod2p(x):
    """
    Custom implementation of math.fmod
    Returns a value between 0 and 2 * math.pi
    """
    rval = math.fmod(x, 2 * math.pi)
    if rval < 0:
        rval += 2 * math.pi
    return rval


def rvel(r, v):
    """
    Convert state vector to mean orbital elements
    """

    # Needed constants
    XJ3 = -2.53881e-6
    XKE = 0.0743669161331734132  # = (G*M)^(1/2)*(er/min)^(3/2) where G
    CK2 = 5.413079e-4
    A3OVK2 = -1 * XJ3 / CK2

    # Position in Earth radii, velocity in km / (Earth radii) / min
    rr2 = r.value / 6378.135
    vv2 = v.value * 60 / 6378.135

    vk = 1.0 / XKE * vv2  # smult
    h = np.cross(rr2, vk)  # cross
    pl = h @ h
    vz = np.array([0.0, 0.0, 1.0])
    n = np.cross(vz, h)
    n = unitv(n)
    rk = np.linalg.norm(rr2)
    rdotk = (rr2 @ vv2) / rk
    rfdotk = np.linalg.norm(h) * XKE / rk
    temp = (rr2 @ n) / rk
    uk = acose(temp)
    if rr2[2] < 0.0:
        uk = 2 * math.pi - uk
    vz = np.cross(vk, h)
    vy = -1 / rk * rr2
    vec = vz + vy
    ek = np.linalg.norm(vec)
    if ek > 1.0:
        print(ek)
        return  # open orbit
    xnodek = math.atan2(n[1], n[0])
    if xnodek < 0.0:
        xnodek += 2 * math.pi
    temp = np.sqrt(h[0] ** 2 + h[1] ** 2)
    xinck = math.atan2(temp, h[2])
    temp = (vec @ n) / ek
    wk = acose(temp)
    if vec[2] < 0:
        wk = fmod2p(2 * math.pi - wk)
    aodp = pl / (1.0 - ek ** 2)
    xn = XKE * aodp ** (-1.5)

    # In the first loop the osculating elements rk, uk, xnodek, xinck, rdotk,
    # and rfdotk are used as anchors to find the corresponding final SGP4
    # mean elements r, u, xnodeo, xincl, rdot, and rfdot.  Several other final
    # mean values based on these are also found: betal, cosio, sinio, theta2,
    # cos2u, sin2u, x3thm1, x7thm1, x1mth2.  In addition, the osculating values
    # initially held by aodp, pl, and xn are replaced by intermediate
    # (not osculating and not mean) values used by SGP4.  The loop converges
    # on the value of pl in about four iterations.

    # seed value for first loop
    xincl = xinck
    u = uk

    for _ in range(0, 99):
        a2 = pl
        betal = math.sqrt(pl / aodp)
        temp1 = CK2 / pl
        temp2 = temp1 / pl
        cosio = math.cos(xincl)
        sinio = math.sin(xincl)
        sin2u = math.sin(2 * u)
        cos2u = math.cos(2 * u)
        theta2 = cosio * cosio
        x3thm1 = 3 * theta2 - 1
        x1mth2 = 1 - theta2
        x7thm1 = 7 * theta2 - 1
        r = (rk - 0.5 * temp1 * x1mth2 * cos2u) / (1 - 1.5 * temp2 * betal * x3thm1)
        u = uk + 0.25 * temp2 * x7thm1 * sin2u
        xnodeo = xnodek - 1.5 * temp2 * cosio * sin2u
        xincl = xinck - 1.5 * temp2 * cosio * sinio * cos2u
        rdot = rdotk + xn * temp1 * x1mth2 * sin2u
        rfdot = rfdotk - xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1)
        pl = (r * rfdot / XKE) ** 2

        # vis-viva equation
        adop = 1 / (2 / r - (rdot ** 2 + rfdot ** 2) / (XKE ** 2))

        xn = XKE * aodp ** (-1.5)
        if math.fabs(a2 - pl) < 1e-13:
            break

    # The next values are calculated from constants and a combination of mean
    # and intermediate quantities from the first loop.  These values all remain
    # fixed and are used in the second loop.

    # preliminary values for the second loops
    ecose = 1 - r / adop
    esine = r * rdot / (XKE * np.sqrt(aodp))  # needed for Kepler's eqn
    elsq = 1 - pl / adop  # intermediate eccentricity squared
    xlcof = 0.125 * A3OVK2 * sinio * (3 + 5 * cosio) / (1 + cosio)
    aycof = 0.25 * A3OVK2 * sinio
    temp1 = esine / (1 + np.sqrt(1 - elsq))
    cosu = math.cos(u)
    sinu = math.sin(u)

    # The second loop normally converges in about six iterations to the final
    # mean value for the eccentricity, eo.  The mean perigee, omegao, is also
    # determined.  Cosepw and sinepw are found to high accuracy and
    # are used to calculate an intermediate value for the eccentric anomaly,
    # temp2.  Temp2 is then used in Kepler's equation to find an intermediate
    # value for the true longitude, capu.

    # seed values for loop
    eo = np.sqrt(elsq)
    omegao = wk
    axn = eo * math.cos(omegao)

    for _ in range(0, 99):
        a2 = eo
        beta = 1 - eo ** 2
        aynl = aycof / (aodp * beta)
        ayn = eo * math.sin(omegao) + aynl
        cosepw = r * cosu / aodp + axn - ayn * temp1
        sinepw = r * sinu / aodp + ayn + axn * temp1
        axn = cosepw * ecose + sinepw * esine
        ayn = sinepw * ecose - cosepw * esine
        omegao = fmod2p(math.atan2(ayn - aynl, axn))
        # use weighted average to tame instability at high eccentricities
        eo = 0.9 * eo + 0.1 * (axn / math.cos(omegao))
        if eo > 0.999:
            eo = 0.999
        if math.fabs(a2 - eo) < 1e-13:
            break

    capu = math.atan2(sinepw, cosepw) - esine  # Kepler's equation
    xll = temp * xlcof * axn

    # xll adjusts the intermediate true longitude
    # capu, to the mean true longitude, xl
    xl = capu - xll

    xmo = fmod2p(xl - omegao)  # mean anomaly

    # The third loop usually converges after three iterations to the
    # mean semi-major axis, a1, which is then used to find the mean motion, xno.
    a0 = aodp
    a1 = a0
    beta2 = np.sqrt(beta)
    temp = 1.5 * CK2 * x3thm1 / (beta * beta2)
    for _ in range(0, 99):
        a2 = a1
        d0 = temp / (a0 ** 2)
        a0 = aodp * (1.0 - d0)
        d1 = temp / (a1 ** 2)
        a1 = a0 / (1 - d1 / 3 - d1 ** 2 - 134 * d1 ** 3 / 81)
        if math.fabs(a2 - a1) < 1e-3:
            break

    xno = XKE * a1 ** (-1.5)

    return xincl, xnodeo, eo, omegao, xmo, xno


def el2rv(inc, raan, ecc, argp, mean_anomaly, mean_motion, epoch):
    """
    Converts mean orbital elements to state vector
    """

    time_tle = epoch.jd - 2433281.5
    sat = Satrec()
    sat.sgp4init(
        WGS84,
        "i",
        0,
        time_tle,
        0.0,
        0.0,
        0.0,
        ecc,
        argp,
        inc,
        mean_anomaly,
        mean_motion,
        raan,
    )

    errorCode, rTEME, vTEME = sat.sgp4(epoch.jd1, epoch.jd2)
    if errorCode != 0:
        raise RuntimeError(SGP4_ERRORS[errorCode])

    pTEME = coord.CartesianRepresentation(rTEME * u.km)
    vTEME = coord.CartesianDifferential(vTEME * u.km / u.s)
    svTEME = TEME(pTEME.with_differentials(vTEME), obstime=epoch)

    svITRS = svTEME.transform_to(coord.ITRS(obstime=epoch))

    orb = Orbit.from_coords(Earth, svITRS)

    return orb.r, orb.v


def rv2el(rr, vv, epoch):
    """
    Computes mean orbital elements from state vector
    """

    epoch_time = Time(epoch, format="datetime", scale="utc")

    # SPG4 k-elements from state vector
    inck, raank, ecck, argpk, mAnomalyk, mMotionk = rvel(rr, vv)

    # SPG4 propagation of k-elements to rr', vv'
    pos, vel = el2rv(inck, raank, ecck, argpk, mAnomalyk, mMotionk, epoch_time)
    inc2, raan2, ecc2, argp2, mAnomaly2, mMotion2 = rvel(
        pos, vel
    )  # SPG4 x-elements from state vectors

    # First correction
    incz = 2 * inck - inc2
    raanz = 2 * raank - raan2
    eccz = 2 * ecck - ecc2
    argpz = 2 * argpk - argp2
    mAnomalyz = 2 * mAnomalyk - mAnomaly2
    mMotionz = 2 * mMotionk - mMotion2

    # second correction is a small adjustment to z-elements
    pos, vel = el2rv(incz, raanz, eccz, argpz, mAnomalyz, mMotionz, epoch_time)
    inc3, raan3, ecc3, argp3, mAnomaly3, mMotion3 = rvel(pos, vel)

    inc = incz + inck - inc3
    raan = raanz + raank - raan3
    ecc = eccz + ecck - ecc3
    argp = argpz + argpk - argp3
    mAnomaly = mAnomalyz + mAnomalyk - mAnomaly3
    mMotion = mMotionz + mMotionk - mMotion3

    return inc, raan, ecc, argp, mAnomaly, mMotion


if __name__ == "__main__":

    # Display some initial data
    print(f" Orbit: {iss}")
    print(" State vector [poliastro]")
    print(f"     r = {iss.r}")
    print(f"     v = {iss.v}")
    print()

    # Reference epoch
    epoch = Time(iss.epoch, format="datetime", scale="utc")

    # Store poliastro orbital elements (osculating)
    ecc_anomaly = angles.nu_to_E(iss.nu, iss.ecc)
    mean_anomaly = angles.E_to_M(ecc_anomaly, iss.ecc)

    # Compute orbital elements required by spg4 (mean)
    inc, raan, ecc, argp, m_ano, m_mot = rv2el(iss.r, iss.v, iss.epoch)

    # Display differences
    print("                Poliastro(osc)              rv2el(mean)")
    print(f"Ecc            :    {iss.ecc:10.5f}{'':15}{ecc:10.5f}")
    print(
        f"Incl  [deg]    :    {math.degrees(iss.inc.value):10.5f}{'':15}{math.degrees(inc):10.5f}"
    )
    print(
        f"n [deg/min]    :    {math.degrees(iss.n.to(u.rad/u.minute).value):10.5f}{'':15}{math.degrees(m_mot):10.5f}"
    )
    print(
        f"RAAN  [deg]    :    {math.degrees(iss.raan.value):10.5f}{'':15}{math.degrees(raan):10.5f}"
    )
    print(
        f"Argp + M [deg] :    {math.degrees(iss.argp.value+mean_anomaly.value):10.5f}{'':15}{math.degrees(argp+m_ano):10.5f}"
    )
    print()

    # Obtain state vector from spg4 and mean elements
    sat = Satrec()
    sat.sgp4init(
        WGS84,
        "i",
        0,
        epoch.jd - 2433281.5,
        0.0,
        0.0,
        0.0,
        ecc,
        argp,
        inc,
        m_ano,
        m_mot,
        raan,
    )
    errorCode, rTEME, vTEME = sat.sgp4(epoch.jd1, epoch.jd2)
    if errorCode != 0:
        raise RuntimeError(SGP4_ERRORS[errorCode])

    # Convert state vector from TEME (True Equator Mean Equinox) to ITRS
    pTEME = coord.CartesianRepresentation(rTEME * u.km)
    vTEME = coord.CartesianDifferential(vTEME * u.km / u.s)
    svTEME = TEME(pTEME.with_differentials(vTEME), obstime=iss.epoch)
    svITRS = svTEME.transform_to(coord.ITRS(obstime=iss.epoch))
    sv = Orbit.from_coords(Earth, svITRS)

    # Display results
    print("State vector [rv2el]")
    print(f"     r = {sv.r}")
    print(f"     v = {sv.v}")
    print()

    print("State vector differences [poliastro - rv2el]")
    print(f"    dr = {iss.r - sv.r}")
    print(f"    dv = {iss.v - sv.v}")
    print()
