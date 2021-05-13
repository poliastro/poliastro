"""
-------------------------------------------------------------------- 
---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ----------
--------------------------------------------------------------------

This file has been ported from the NRLMSISE-00 C source code package 
- release 20041227

The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
Doug Drob. Model is also available as a NRLMSISE-00 distribution 
package in FORTRAN (link to model FORTRAN: 
https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/)

Dominik Brodowski implemented and maintains this C version. You can
reach him at mail@brodo.de. 

Version Dated: 2019-07-09 1255 hrs

This is the Main Code for NRLMSISE00 Model with class implementations for input,
flags & output.
"""

import math

import astropy
import numpy as np
from astropy import units as u

from poliastro.earth.atmosphere.nrlmsise00_data import (
    pavgm,
    pd,
    pdl,
    pdm,
    pma,
    ps,
    pt,
    ptl,
    ptm,
    sam,
)

# Constant Values (Units to be added later)
dgtr = 1.74533e-2
rgas = 831.4
sr = 7.2722e-5
dr = 1.72142e-2
hr = 0.2618
pset = 2.0


class nrlmsise_flags(object):
    """
    Switches: to turn on and off particular variations use these switches.
    0 is off, 1 is on, and 2 is main effects off but cross terms on.

    Standard values are 0 for switch 0 and 1 for switches 1 to 23. The
    array "switches" needs to be set accordingly by the calling program.
    The arrays sw and swc are set internally.

    Parameters
    ----------
    switches: list
        i , explanation , size of list = 24
    0  : int (0,1,2)
        Output in meters and kilograms (SI Units) instead of centimeters and grams (CGS)
    1  : int (0,1,2)
        F10.7 Effect on Mean
    2  : int (0,1,2)
        Time independent
    3  : int (0,1,2)
        Symmetrical Annual
    4  : int (0,1,2)
        Symmetrical Semiannual
    5  : int (0,1,2)
        Asymmetrical Annual
    6  : int (0,1,2)
        asymmetrical semiannual
    7  : int (0,1,2)
        Diurnal
    8  : int (0,1,2)
        Semidiurnal
    9  : int (0,1,2)
        Daily AP [when this is set to -1 (!) the ap_a in nrlmsise_input must
                    point to an instance of class ap_array]
    10 : int (0,1,2)
        All UT/Long Effects
    11 : int (0,1,2)
        Longitudinal
    12 : int (0,1,2)
        UT and Mixed UT/long
    13 : int (0,1,2)
        Mixed AP/UT/LONG
    14 : int (0,1,2)
        Terdiurnal
    15 : int (0,1,2)
        Departures from Diffusive Equilibrium
    16 : int (0,1,2)
        All TINF var
    17 : int (0,1,2)
        All TLB var
    18 : int (0,1,2)
        All TN1 var
    19 : int (0,1,2)
        All S var
    20 : int (0,1,2)
        All TN2 var
    21 : int (0,1,2)
        All NLB var
    22 : int (0,1,2)
        All TN3 var
    23 : int (0,1,2)
        Turbo Scale Height var
    """

    def __init__(self):
        self.switches = [0 for _ in range(24)]
        self.sw = [0.0 for _ in range(24)]
        self.swc = [0.0 for _ in range(24)]


class ap_array:
    """
    Array containing the following magnetic values:

    Parameters
    ----------
    0 : float
        Daily AP
    1 : float
        3 hr AP index for current time
    2 : float
        3 hr AP index for 3 hrs before current time
    3 : float
        3 hr AP index for 6 hrs before current time
    4 : float
        3 hr AP index for 9 hrs before current time
    5 : float
        Average of eight 3 hr AP indicies from 12 to 33 hrs
        prior to current time
    6 : float
        Average of eight 3 hr AP indicies from 36 to 57 hrs
        prior to current time
    """

    def __init__(self):
        self.a = [0.0 for _ in range(7)]


class nrlmsise_input:
    """
    NOTES ON INPUT VARIABLES:

    UT, Local Time, and Longitude are used independently in the
    model and are not of equal importance for every situation.
    For the most physically realistic calculation these three
    variables should be consistent (lst=sec/3600 + g_long/15).
    The Equation of Time departures from the above formula
    for apparent local time can be included if available but
    are of minor importance.

    f107 and f107A values used to generate the model correspond
    to the 10.7 cm radio flux at the actual distance of the Earth
    from the Sun rather than the radio flux at 1 AU. The following
    site provides both classes of values:
    ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/

    f107, f107A, and ap effects are neither large nor well
    established below 80 km and these parameters should be set to
    150., 150., and 4. respectively.

    Parameters
    ----------
    year : int
        Year for computation
    doy : int
        Day of the Year
    sec : float
        Seconds of the day (UT)
    alt : float
        Altitude (km)
    g_lat : float
        Geodetic Latitude
    g_long : float
        Geodetic Longitude
    lst : float
        Local Apparent Solar Time (hrs)
    F107A : float
        81 Day Average for 10.7 Flux (Centered on day)
        Defaults to 150.0
    F107 : float
        Daily F10.7 Flux for previous day
        Defaults to 150.0
    ap : float
        Magnetic Index (Daily)
        Defaults to 4.0
    ap_a : ~poliastro.earth.atmosphere.nrlmsise00.ap_array
        Check `ap_array` for information
    """

    def __init__(
        self,
        year=0,
        doy=0,
        sec=0.0,
        alt=0.0,
        g_lat=0.0,
        g_long=0.0,
        lst=0.0,
        f107A=150.0,
        f107=150.0,
        ap=4.0,
        ap_a=None,
    ):
        self.year = year
        self.doy = doy
        self.sec = sec
        self.alt = alt
        self.g_lat = g_lat
        self.g_long = g_long
        self.lst = lst
        self.f107A = f107A
        self.f107 = f107
        self.ap = ap
        self.ap_array = ap_a

    def set_lst(self):
        """
        Generates the value for Local Apparent Solar Time using the formula
        lst = sec/3600 + g_long/15

        Parameters
        ----------
        sec : float
            Seconds of the Day
        g_long : float
            Geodetic Longitude

        Returns
        -------
        lst : float
            Local Apparent Solar Time (UT) calculated using longitude value
        """
        self.lst = (self.sec / 3600) + (self.g_long / 15)
        return self.lst


class nrlmsise_output:
    """
    Output Class for storing values for NRLMSISE00 Computation

    O, H, and N are set to zero below 72.5 km

    t[0], Exospheric temperature, is set to global average for
    altitudes below 120 km. The 120 km gradient is left at global
    average value for altitudes below 72 km.

    d[5], Total Mass Density, is NOT the same for subroutines GTD7
    and GTD7D

    SUBROUTINE GTD7 : d[5] is the sum of the mass densities of the
    species labeled by indices 0-4 and 6-7 in output variable d.
    This includes He, O, N2, O2, Ar, H, and N but does NOT include
    anomalous oxygen (species index 8).

    SUBROUTINE GTD7D : d[5] is the "effective total mass density
    for drag" and is the sum of the mass densities of all species
    in this model, INCLUDING anomalous oxygen.

    Parameters
    ----------
    d[0] : float
        HE Number Density(cm-3)
    d[1] : float
        O Number Density(cm-3)
    d[2] : float
        N2 Number Density(cm-3)
    d[3] : float
        O2 Number Density(cm-3)
    d[4] : float
        AR Number Density(cm-3)
    d[5] : float
        Total Mass denisty(g/cm3) [includes d[8] in td7d]
    d[6] : float
        H Number Density(cm-3)
    d[7] : float
        N Number Density(cm-3)
    d[8] : float
        Anomalous Oxygen Number Density(cm-3)
    t[0] : float
        Exospheric Temperature
    t[1] : float
        Temperature at altitude
    """

    def __init__(self):
        self.d = [0.0 for _ in range(9)]  # Densities
        self.t = [0.0 for _ in range(2)]  # Temperatures


# Values to set (global variables)
# PARMB
gsurf = [0.0]
re = [0.0]

# GTS3C
dd = 0.0

# DMIX
dm04 = 0.0
dm16 = 0.0
dm28 = 0.0
dm32 = 0.0
dm40 = 0.0
dm01 = 0.0
dm14 = 0.0

# MESO7
meso_tn1 = [0.0 for _ in range(5)]
meso_tn2 = [0.0 for _ in range(4)]
meso_tn3 = [0.0 for _ in range(5)]
meso_tgn1 = [0.0 for _ in range(2)]
meso_tgn2 = [0.0 for _ in range(2)]
meso_tgn3 = [0.0 for _ in range(2)]

# LPOLY
dfa = 0.0
plg = [[0.0 for _ in range(9)] for _ in range(4)]
ctloc = 0.0
stloc = 0.0
c2tloc = 0.0
s2tloc = 0.0
s3tloc = 0.0
c3tloc = 0.0
apdf = 0.0
apt = [0.0 for _ in range(4)]

# TSELEC
def tselec(flags):
    """
    Parameters
    ----------
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Flags for various values (0,1,2)
        Check class nrlmsise_flags for complete description
    """
    for index in range(0, 24):
        if not index == 9:
            if flags.switches[index] == 1:
                flags.sw[index] = 1
            else:
                flags.sw[index] = 0

            if flags.switches[index] > 0:
                flags.swc[index] = 1
            else:
                flags.swc[index] = 0
        else:
            flags.sw[index] = flags.switches[index]
            flags.swc[index] = flags.switches[index]
    return


# GLATF
def glatf(lat, gv, reff):
    """
    Parameters
    ----------
    lat : float
        Latitude Value
    gv : float
        Effective Gravity
    reff : float
        Effective Radius of Earth
    """
    c2 = np.cos(2.0 * dgtr * lat)
    gv[0] = 980.616 * (1.0 - 0.0026373 * c2)
    reff[0] = 2.0 * gv[0] / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5
    return


def ccor(alt, r, h1, zh):
    """
    Chemistry/Dissociation Correction for MSIS Models

    Parameters
    ----------
    alt : float
        Altitude
    r   : float
        Target Ratio
    h1  : float
        Transition Scale Length
    zh  : float
        Altitude of (1/2 * r)

    Returns
    -------
    exp value : float
        Exponential value based on Altitude
    """
    e = (alt - zh) / h1
    if e > 70:
        return np.exp(0)
    if e < -70:
        return np.exp(r)
    ex = np.exp(e)
    e = r / (1.0 + ex)
    return np.exp(e)


def ccor2(alt, r, h1, zh, h2):
    """
    Chemistry/Dissociation Correction for MSIS Models

    alt : float
        Altitude
    r   : float
        Target Ratio
    h1  : float
        Transition Scale Length
    zh  : float
        Altitude of (1/2 * r)
    h2  : float
        Transition Scale Length 2

    Returns
    -------
    exp value : float
        Exponential value based on Altitude
    """
    e1 = (alt - zh) / h1
    e2 = (alt - zh) / h2
    if (e1 > 70) or (e2 > 70):
        return np.exp(0)
    if (e1 < -70) and (e2 < -70):
        return np.exp(r)
    ex1 = np.exp(e1)
    ex2 = np.exp(e2)
    ccor2v = r / (1.0 + 0.5 * (ex1 + ex2))
    return np.exp(ccor2v)


def scalh(alt, xm, temp):
    """
    Parameters
    ----------
    alt : float
        Altitude
    xm : float
        Species Molecular Weight
    temp : float
        Temperature
    """
    g = gsurf[0] / math.pow((1.0 + alt / re[0]), 2.0)
    g = rgas * temp / (g * xm)
    return g


def dnet(dd, dm, zhm, xmm, xm):
    """
    TurboPause Correction For MSIS Model
    Root Mean Density

    Parameters
    ----------
    dd : float
        Diffusive Density
    dm : float
        Full Mixed Density
    zhm : float
        Transition Scale Length
    xmm : float
        Full Mixed Molecular Weight
    xm : float
        Species Molecular Weight
    dnet : float
        Combined Density

    Returns
    -------
    dm : float
        if (dd = 0) , Full Mixed Density
    dd : float
        if (dd = 0) , Diffusive Density
    a : float
    """
    a = zhm / (xmm - xm)
    if not ((dm > 0) and (dd > 0)):
        # print(f"dnet log error {dm}, {dd}, {xm} \n")
        if (dd == 0) and (dm == 0):
            dd = 1
        if dm == 0:
            return dd
        if dd == 0:
            return dm

    ylog = a * math.log(dm / dd)
    if ylog < -10:
        return dd
    if ylog > 10:
        return dm
    a = dd * math.pow((1.0 + np.exp(ylog)), (1.0 / a))
    return a


def splini(xa, ya, y2a, n, x, y):
    """
    Integrate Cubic Spline Function from XA(1) to X

    Parameters
    ----------
    xa : list
        Array of Tabulated Functions in Ascending Order by x
    ya  : list
        Array of Tabulated Functions in Ascending Order by x
    y2a : list
        Array of Second Derivatives
    n : int
        Size of Arrays xa, ya, y2a
    x : float
        Abscissa Endpoint for Integration
    y : float
        Output Value
    """
    yi = 0
    klo = 0
    khi = 1
    while (x > xa[klo]) and (khi < n):
        xx = x
        if khi < (n - 1):
            if x < xa[khi]:
                xx = x
            else:
                xx = xa[khi]
        h = xa[khi] - xa[klo]
        a = (xa[khi] - xx) / h
        b = (xx - xa[klo]) / h
        a2 = a * a
        b2 = b * b
        yi += (
            (1.0 - a2) * ya[klo] / 2.0
            + b2 * ya[khi] / 2.0
            + (
                (-(1.0 + a2 * a2) / 4.0 + a2 / 2.0) * y2a[klo]
                + (b2 * b2 / 4.0 - b2 / 2.0) * y2a[khi]
            )
            * h
            * h
            / 6.0
        ) * h
        klo += 1
        khi += 1
    y[0] = yi
    return


def splint(xa, ya, y2a, n, x, y):
    """
    Calculate Cubic Spline Interp Value
    Adapted From Numerical Recipes by PRESS ET AL.

    Parameters
    ----------
    xa : list
        Array of Tabulated Functions in Ascending Order by x
    ya  : list
        Array of Tabulated Functions in Ascending Order by x
    y2a : list
        Array of Second Derivatives
    n : int
        Size of Arrays xa, ya, y2a
    x : float
        Abscissa for Interpolation
    y : float
        Output Value
    """
    klo = 0
    khi = n - 1
    while (khi - klo) > 1:
        k = int((khi + klo) / 2)
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    if h == 0.0:
        print("bad XA input to Splint")
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    yi = (
        a * ya[klo]
        + b * ya[khi]
        + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0
    )
    y[0] = yi
    return


def spline(x, y, n, yp1, ypn, y2):
    """
    Calculate 2nd Derivatives of Cubic Spline Interp Function
    Adapted from Numerical Recipes by PRESS ET AL.

    Parameters
    ----------
    x : list
        Arrays of Tabulated Function in Ascending Order by x
    y : list
        Arrays of Tabulated Function in Ascending Order by x
    n : int
        Sizes of Arrays x,y
    yp1 : float
        Specified Derivative value at x[0]
    ypn : float
        Specified Derivative value at x[n-1]
        Values >= 1E30 Signal Signal Second Derivative Zero
    y2 : list
        Output Array of Second Derivatives
    """
    u = [0.0 for _ in range(n)]
    if yp1 > 0.99e30:
        y2[0] = 0
        u[0] = 0
    else:
        y2[0] = -0.5
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1)
    for index in range(1, n - 1):
        sig = (x[index] - x[index - 1]) / (x[index + 1] - x[index - 1])
        p = sig * y2[index - 1] + 2.0
        y2[index] = (sig - 1.0) / p
        u[index] = (
            6.0
            * (
                (y[index + 1] - y[index]) / (x[index + 1] - x[index])
                - (y[index] - y[index - 1]) / (x[index] - x[index - 1])
            )
            / (x[index + 1] - x[index - 1])
            - sig * u[index - 1]
        ) / p
    if ypn > 0.99e30:
        qn = 0
        un = 0
    else:
        qn = 0.5
        un = (3.0 / (x[n - 1] - x[n - 2])) * (
            ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
        )
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)
    index = n - 2
    while index >= 0:
        y2[index] = y2[index] * y2[index + 1] + u[index]
        index -= 1
    return


def zeta(zz, zl):
    """
    Parameters
    ----------
    zz : float
    zl : float

    Returns
    -------
    zeta_value : float
        Geopotential Difference between the parameters
        Computation based on the following formula
        (zz - zl) * (re + zl) / (re + zz)
    """
    return (zz - zl) * (re[0] + zl) / (re[0] + zz)


def densm(alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2):
    """
    Calculate Temperature & Density Profiles for Lower Atoms.

    Parameters
    ----------
    alt : float
        Altitude
    d0 : float
    xm : float
        Species Molecular Weight
    tz : list
        Contains only one value
    mn3 : list
        Troposphere / Stratosphere
    zn3 : list
        Troposphere / Stratosphere
    tn3 : list
        Troposphere / Stratosphere
    tgn3 : list
        Troposphere / Stratosphere
    mn2 : list
        Stratosphere / Mesosphere
    zn2 : list
        Stratosphere / Mesosphere
    tn2 : list
        Stratosphere / Mesosphere
    tgn2 : list
        Stratosphere / Mesosphere

    Returns
    -------
    tz : float
    densm_temp : float
    """
    xs = [0.0 for _ in range(10)]
    ys = [0.0 for _ in range(10)]
    y2out = [0.0 for _ in range(10)]
    densm_temp = d0
    if alt > zn2[0]:
        if xm == 0.0:
            return tz[0]
        else:
            return d0

    # Stratosphere / Mesosphere Temperature
    if alt > zn2[mn2 - 1]:
        z = alt
    else:
        z = zn2[mn2 - 1]
    mn = mn2
    z1 = zn2[0]
    z2 = zn2[mn - 1]
    t1 = tn2[0]
    t2 = tn2[mn - 1]
    zg = zeta(z, z1)
    zgdif = zeta(z2, z1)

    # Set Up Spline Nodes
    for i in range(0, mn):
        xs[i] = zeta(zn2[i], z1) / zgdif
        ys[i] = 1.0 / tn2[i]
    yd1 = -tgn2[0] / (t1 * t1) * zgdif
    yd2 = -tgn2[1] / (t2 * t2) * zgdif * pow(((re[0] + z2) / (re[0] + z1)), 2.0)

    # Calculate Spline Coefficients
    spline(xs, ys, mn, yd1, yd2, y2out)
    x = zg / zgdif
    y = [0.0]
    splint(xs, ys, y2out, mn, x, y)

    # Temperature at Altitude
    tz[0] = 1.0 / y[0]
    if not xm == 0.0:
        # Calculates Stratosphere / Mesosphere Density
        glb = gsurf[0] / pow((1.0 + z1 / re[0]), 2.0)
        gamm = xm * glb * zgdif / rgas

        # Integrate Temperate Profile
        yi = [0.0]
        splini(xs, ys, y2out, mn, x, yi)
        expl = gamm * yi[0]
        if expl > 50.0:
            expl = 50.0

        # Density at Altitude
        densm_temp = densm_temp * (t1 / tz[0]) * np.exp(-expl)

    if alt > zn3[0]:
        if xm == 0.0:
            return tz[0]
        else:
            return densm_temp

    # Troposphere / Stratosphere Temperature
    z = alt
    mn = mn3
    z1 = zn3[0]
    z2 = zn3[mn - 1]
    t1 = tn3[0]
    t2 = tn3[mn - 1]
    zg = zeta(z, z1)
    zgdif = zeta(z2, z1)

    # Set Up Spline Nodes
    for i in range(0, mn):
        xs[i] = zeta(zn3[i], z1) / zgdif
        ys[i] = 1.0 / tn3[i]
    yd1 = -tgn3[0] / (t1 * t1) * zgdif
    yd2 = -tgn3[1] / (t2 * t2) * zgdif * pow(((re[0] + z2) / (re[0] + z1)), 2.0)

    # Calculate Spline Coefficients
    spline(xs, ys, mn, yd1, yd2, y2out)
    x = zg / zgdif
    y = [0.0]
    splint(xs, ys, y2out, mn, x, y)

    # Temperature at Altitude
    tz[0] = 1.0 / y[0]
    if not xm == 0.0:
        # Calculate Tropospheric / Stratospheric Density
        glb = gsurf[0] / pow((1.0 + z1 / re[0]), 2.0)
        gamm = xm * glb * zgdif / rgas

        # Integrate Temperature Profile
        yi = [0.0]
        splini(xs, ys, y2out, mn, x, yi)
        expl = gamm * yi[0]
        if expl > 50.0:
            expl = 50.0

        # Density at Altitude
        densm_temp = densm_temp * (t1 / tz[0]) * np.exp(-expl)

    if xm == 0.0:
        return tz[0]
    else:
        return densm_temp


def densu(alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, mn1, zn1, tn1, tgn1):
    """
    Calculate Temperature & Density Profile for MSIS Model

    Parameters
    ----------
    alt : float
        Altitude
    dlb : float
    tinf : float
        Temperature at inf
    tlb : float
        Temperature at
    xm : float
        Species Molecular Weight
    alpha : float
    tz : list
        Contains only one value
    zlb : float
    s2 : float
    mn1 : list
    zn1 : list
    tn1 : list
    tgn1 : list

    Returns
    -------
    densu_temp : float
    """
    # New Lower Thermo Polynomial
    densu_temp = 1.0
    xs = [0.0 for _ in range(5)]
    ys = [0.0 for _ in range(5)]
    y2out = [0.0 for _ in range(5)]

    # Joining Altitudes of Bates & Spline
    za = zn1[0]
    if alt > za:
        z = alt
    else:
        z = za

    # Geopotential Altitude Difference from ZBL
    zg2 = zeta(z, zlb)

    # Bates Temperature
    tt = tinf - (tinf - tlb) * np.exp(-s2 * zg2)
    ta = tt
    tz[0] = tt
    densu_temp = tz[0]

    if alt < za:
        # Calculates Temperature Below ZA
        # Temperature Gradient at ZA from Bates Profile.
        dta = (tinf - ta) * s2 * pow(((re[0] + zlb) / (re[0] + za)), 2.0)
        tgn1[0] = dta
        tn1[0] = ta
        if alt > zn1[mn1 - 1]:
            z = alt
        else:
            z = zn1[mn1 - 1]
        mn = mn1
        z1 = zn1[0]
        z2 = zn1[mn - 1]
        t1 = tn1[0]
        t2 = tn1[mn - 1]

        # Geopotential Difference from Z1
        zg = zeta(z, z1)
        zgdif = zeta(z2, z1)

        # Set Up Spline Nodes
        for i in range(0, mn):
            xs[i] = zeta(zn1[i], z1) / zgdif
            ys[i] = 1.0 / tn1[i]

        # End Node Derivatives
        yd1 = -tgn1[0] / (t1 * t1) * zgdif
        yd2 = -tgn1[1] / (t2 * t2) * zgdif * pow(((re[0] + z2) / (re[0] + z1)), 2.0)

        # Calculates Spline Coefficients
        spline(xs, ys, mn, yd1, yd2, y2out)
        x = zg / zgdif
        y = [0.0]
        splint(xs, ys, y2out, mn, x, y)

        # Temperature at Altitude
        tz[0] = 1.0 / y[0]
        densu_temp = tz[0]

    if xm == 0:
        return densu_temp

    # Calculate Density Above ZA
    glb = gsurf[0] / pow((1.0 + zlb / re[0]), 2.0)
    gamma = xm * glb / (s2 * rgas * tinf)
    expl = np.exp(-s2 * gamma * zg2)
    if expl > 50.0:
        expl = 50.0
    if tt <= 0:
        expl = 50.0

    # Denstiy at Altitude
    densa = dlb * pow((tlb / tt), (1.0 + alpha + gamma)) * expl
    densu_temp = densa
    if alt >= za:
        return densu_temp

    # Calculate Density Below ZA
    glb = gsurf[0] / pow((1.0 + z1 / re[0]), 2.0)
    gamm = xm * glb * zgdif / rgas

    # Integrate Spline Temperatures
    yi = [0]
    splini(xs, ys, y2out, mn, x, yi)
    expl = gamm * yi[0]
    if expl > 50.0:
        expl = 50.0
    if tz[0] <= 0:
        expl = 50.0

    # Density at Altitude
    densu_temp = densu_temp * pow((t1 / tz[0]), (1.0 + alpha)) * np.exp(-expl)
    return densu_temp


# GLOBE7

# 3 hr Magnetic Activity Function
# Eq. A24d
def g0(a, p):
    """
    Parameters
    ----------
    a : float
    p : float

    Returns
    -------
    """
    return (
        a
        - 4.0
        + (p[25] - 1.0)
        * (
            a
            - 4.0
            + (np.exp(-math.sqrt(p[24] * p[24]) * (a - 4.0)) - 1.0)
            / math.sqrt(p[24] * p[24])
        )
    )


# Eq. A24c
def sumex(ex):
    """
    Parameters
    ----------
    ex : float

    Returns
    -------
    """
    return 1.0 + (1.0 - pow(ex, 19.0)) / (1.0 - ex) * pow(ex, 0.5)


# Eq. A24a
def sg0(ex, p, ap):
    """
    Parameters
    ----------
    ex ([type]): [description]
    p ([type]): [description]
    ap ([type]): [description]

    Returns
    -------
    """
    return (
        g0(ap[1], p)
        + (
            g0(ap[2], p) * ex
            + g0(ap[3], p) * ex * ex
            + g0(ap[4], p) * pow(ex, 3.0)
            + (g0(ap[5], p) * pow(ex, 4.0) + g0(ap[6], p) * pow(ex, 12.0))
            * (1.0 - pow(ex, 8.0))
            / (1.0 - ex)
        )
    ) / sumex(ex)


def globe7(p, input, flags):
    """
    Calculate G(L) Function

    Parameters
    ----------
    p : list
    input : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_input
        Nrlmsise_input class instance
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Nrlmsise_flags class instance

    Returns
    -------
    tinf : float
    """
    # Upper Thermosphere Parameters
    t = [0 for _ in range(15)]
    sw9 = 1
    tloc = input.lst

    if flags.sw[9] > 0:
        sw9 = 1
    elif flags.sw[9] < 0:
        sw9 = -1
    xlong = input.g_long

    # Calculate Legendre Polynomials
    c = np.sin(input.g_lat * dgtr)
    s = np.cos(input.g_lat * dgtr)
    c2 = c * c
    c4 = c2 * c2
    s2 = s * s

    plg[0][1] = c
    plg[0][2] = 0.5 * (3.0 * c2 - 1.0)
    plg[0][3] = 0.5 * (5.0 * c * c2 - 3.0 * c)
    plg[0][4] = (35.0 * c4 - 30.0 * c2 + 3.0) / 8.0
    plg[0][5] = (63.0 * c2 * c2 * c - 70.0 * c2 * c + 15.0 * c) / 8.0
    plg[0][6] = (11.0 * c * plg[0][5] - 5.0 * plg[0][4]) / 6.0

    plg[1][1] = s
    plg[1][2] = 3.0 * c * s
    plg[1][3] = 1.5 * (5.0 * c2 - 1.0) * s
    plg[1][4] = 2.5 * (7.0 * c2 * c - 3.0 * c) * s
    plg[1][5] = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s
    plg[1][6] = (11.0 * c * plg[1][5] - 6.0 * plg[1][4]) / 5.0

    plg[2][2] = 3.0 * s2
    plg[2][3] = 15.0 * s2 * c
    plg[2][4] = 7.5 * (7.0 * c2 - 1.0) * s2
    plg[2][5] = 3.0 * c * plg[2][4] - 2.0 * plg[2][3]
    plg[2][6] = (11.0 * c * plg[2][5] - 7.0 * plg[2][4]) / 4.0
    plg[2][7] = (13.0 * c * plg[2][6] - 8.0 * plg[2][5]) / 5.0
    plg[3][3] = 15.0 * s2 * s
    plg[3][4] = 105.0 * s2 * s * c
    plg[3][5] = (9.0 * c * plg[3][4] - 7.0 * plg[3][3]) / 2.0
    plg[3][6] = (11.0 * c * plg[3][5] - 8.0 * plg[3][4]) / 3.0

    if not (((flags.sw[7] == 0) and (flags.sw[8] == 0)) and (flags.sw[14] == 0)):
        global stloc
        stloc = np.sin(hr * tloc)
        global ctloc
        ctloc = np.cos(hr * tloc)
        global s2tloc
        s2tloc = np.sin(2.0 * hr * tloc)
        global c2tloc
        c2tloc = np.cos(2.0 * hr * tloc)
        global s3tloc
        s3tloc = np.sin(3.0 * hr * tloc)
        global c3tloc
        c3tloc = np.cos(3.0 * hr * tloc)

    cd32 = np.cos(dr * (input.doy - p[31]))
    cd18 = np.cos(2.0 * dr * (input.doy - p[17]))
    cd14 = np.cos(dr * (input.doy - p[13]))
    cd39 = np.cos(2.0 * dr * (input.doy - p[38]))
    p32 = p[31]
    p18 = p[17]
    p14 = p[13]
    p39 = p[38]

    # F10.7 EFFECT
    df = input.f107 - input.f107A
    global dfa
    dfa = input.f107A - 150.0
    t[0] = (
        p[19] * df * (1.0 + p[59] * dfa)
        + p[20] * df * df
        + p[21] * dfa
        + p[29] * pow(dfa, 2.0)
    )
    f1 = 1.0 + (p[47] * dfa + p[19] * df + p[20] * df * df) * flags.swc[1]
    f2 = 1.0 + (p[49] * dfa + p[19] * df + p[20] * df * df) * flags.swc[1]

    # Time Independent
    t[1] = (
        (p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6])
        + (p[14] * plg[0][2]) * dfa * flags.swc[1]
        + p[26] * plg[0][1]
    )

    # Symmetrical Annual
    t[2] = p[18] * cd32

    # Symmetrical Semi-Annual
    t[3] = (p[15] + p[16] * plg[0][2]) * cd18

    # Asymmetrical Annual
    t[4] = f1 * (p[9] * plg[0][1] + p[10] * plg[0][3]) * cd14

    # Asymmetrical Semi-Annual
    t[5] = p[37] * plg[0][1] * cd39

    # DIURNAL
    if flags.sw[7]:
        t71 = (p[11] * plg[1][2]) * cd14 * flags.swc[5]
        t72 = (p[12] * plg[1][2]) * cd14 * flags.swc[5]
        t[6] = f2 * (
            (p[3] * plg[1][1] + p[4] * plg[1][3] + p[27] * plg[1][5] + t71) * ctloc
            + (p[6] * plg[1][1] + p[7] * plg[1][3] + p[28] * plg[1][5] + t72) * stloc
        )

    # Semi-DIURNAL
    if flags.sw[8]:
        t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * flags.swc[5]
        t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * flags.swc[5]
        t[7] = f2 * (
            (p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc
            + (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc
        )

    # TerDIURNAL
    if flags.sw[14]:
        t[13] = f2 * (
            (
                p[39] * plg[3][3]
                + (p[93] * plg[3][4] + p[46] * plg[3][6]) * cd14 * flags.swc[5]
            )
            * s3tloc
            + (
                p[40] * plg[3][3]
                + (p[94] * plg[3][4] + p[48] * plg[3][6]) * cd14 * flags.swc[5]
            )
            * c3tloc
        )

    # Magnetic activity based on daily ap
    if flags.sw[9] == -1:
        ap = input.ap_a
        if p[51] != 0:
            exp1 = np.exp(
                -10800.0
                * math.sqrt(p[51] * p[51])
                / (1.0 + p[138] * (45.0 - math.sqrt(input.g_lat * input.g_lat)))
            )
            if exp1 > 0.99999:
                exp1 = 0.99999
            if p[24] < 1.0e-4:
                p[24] = 1.0e-4
            apt[0] = sg0(exp1, p, ap.a)
            # /* apt[1]=sg2(exp1,p,ap->a);
            #   apt[2]=sg0(exp2,p,ap->a);
            #   apt[3]=sg2(exp2,p,ap->a);
            # */
            if flags.sw[9]:
                t[8] = apt[0] * (
                    p[50]
                    + p[96] * plg[0][2]
                    + p[54] * plg[0][4]
                    + (p[125] * plg[0][1] + p[126] * plg[0][3] + p[127] * plg[0][5])
                    * cd14
                    * flags.swc[5]
                    + (p[128] * plg[1][1] + p[129] * plg[1][3] + p[130] * plg[1][5])
                    * flags.swc[7]
                    * np.cos(hr * (tloc - p[131]))
                )

    else:
        apd = input.ap - 4.0
        p44 = p[43]
        p45 = p[44]
        if p44 < 0:
            p44 = 1.0e-5
        global apdf
        apdf = apd + (p45 - 1.0) * (apd + (np.exp(-p44 * apd) - 1.0) / p44)
        if flags.sw[9]:
            t[8] = apdf * (
                p[32]
                + p[45] * plg[0][2]
                + p[34] * plg[0][4]
                + (p[100] * plg[0][1] + p[101] * plg[0][3] + p[102] * plg[0][5])
                * cd14
                * flags.swc[5]
                + (p[121] * plg[1][1] + p[122] * plg[1][3] + p[123] * plg[1][5])
                * flags.swc[7]
                * np.cos(hr * (tloc - p[124]))
            )

    if (flags.sw[10]) and (input.g_long > -1000.0):
        # Longitudinal
        if flags.sw[11]:
            t[10] = (1.0 + p[80] * dfa * flags.swc[1]) * (
                (
                    p[64] * plg[1][2]
                    + p[65] * plg[1][4]
                    + p[66] * plg[1][6]
                    + p[103] * plg[1][1]
                    + p[104] * plg[1][3]
                    + p[105] * plg[1][5]
                    + flags.swc[5]
                    * (p[109] * plg[1][1] + p[110] * plg[1][3] + p[111] * plg[1][5])
                    * cd14
                )
                * np.cos(dgtr * input.g_long)
                + (
                    p[90] * plg[1][2]
                    + p[91] * plg[1][4]
                    + p[92] * plg[1][6]
                    + p[106] * plg[1][1]
                    + p[107] * plg[1][3]
                    + p[108] * plg[1][5]
                    + flags.swc[5]
                    * (p[112] * plg[1][1] + p[113] * plg[1][3] + p[114] * plg[1][5])
                    * cd14
                )
                * np.sin(dgtr * input.g_long)
            )

        # UT and Mixed UT, Longitude
        if flags.sw[12]:
            t[11] = (
                (1.0 + p[95] * plg[0][1])
                * (1.0 + p[81] * dfa * flags.swc[1])
                * (1.0 + p[119] * plg[0][1] * flags.swc[5] * cd14)
                * (
                    (p[68] * plg[0][1] + p[69] * plg[0][3] + p[70] * plg[0][5])
                    * np.cos(sr * (input.sec - p[71]))
                )
            )

            t[11] += (
                flags.swc[11]
                * (p[76] * plg[2][3] + p[77] * plg[2][5] + p[78] * plg[2][7])
                * np.cos(sr * (input.sec - p[79]) + 2.0 * dgtr * input.g_long)
                * (1.0 + p[137] * dfa * flags.swc[1])
            )

        # UT, Longitude Magnetic Activity
        if flags.sw[13]:
            if flags.sw[9] == -1:
                if p[51]:
                    t[12] = (
                        apt[0]
                        * flags.swc[11]
                        * (1.0 + p[132] * plg[0][1])
                        * (
                            (p[52] * plg[1][2] + p[98] * plg[1][4] + p[67] * plg[1][6])
                            * np.cos(dgtr * (input.g_long - p[97]))
                        )
                        + apt[0]
                        * flags.swc[11]
                        * flags.swc[5]
                        * (p[133] * plg[1][1] + p[134] * plg[1][3] + p[135] * plg[1][5])
                        * cd14
                        * np.cos(dgtr * (input.g_long - p[136]))
                        + apt[0]
                        * flags.swc[12]
                        * (p[55] * plg[0][1] + p[56] * plg[0][3] + p[57] * plg[0][5])
                        * np.cos(sr * (input.sec - p[58]))
                    )
            else:
                t[12] = (
                    apdf
                    * flags.swc[11]
                    * (1.0 + p[120] * plg[0][1])
                    * (
                        (p[60] * plg[1][2] + p[61] * plg[1][4] + p[62] * plg[1][6])
                        * np.cos(dgtr * (input.g_long - p[63]))
                    )
                    + apdf
                    * flags.swc[11]
                    * flags.swc[5]
                    * (p[115] * plg[1][1] + p[116] * plg[1][3] + p[117] * plg[1][5])
                    * cd14
                    * np.cos(dgtr * (input.g_long - p[118]))
                    + apdf
                    * flags.swc[12]
                    * (p[83] * plg[0][1] + p[84] * plg[0][3] + p[85] * plg[0][5])
                    * np.cos(sr * (input.sec - p[75]))
                )

    # Parms not used: 82, 89, 99, 139-149
    tinf = p[30]
    for ind in range(14):
        tinf = tinf + abs(flags.sw[ind + 1]) * t[ind]
    return tinf


# GLOB7S
def glob7s(p, input, flags):
    """
    Version of GLOBE for Lower Atmosphere 10/26/99

    Parameters
    ----------
    p : list
    input : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_input
        Nrlmsise_input class instance
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Nrlmsise_flags class instance

    Returns
    -------
    tt : float
    """
    t = [0.0 for _ in range(14)]
    pset = 2.0

    # Confirm Parameter Set
    if p[99] == 0:
        p[99] = pset
    if not p[99] == pset:
        print("Wrong parameter set for GLOB7S")
        return -1

    # for ind in range(0, 14):
    #     t[ind] == 0.0

    cd32 = np.cos(dr * (input.doy - p[31]))
    cd18 = np.cos(2.0 * dr * (input.doy - p[17]))
    cd14 = np.cos(dr * (input.doy - p[13]))
    cd39 = np.cos(2.0 * dr * (input.doy - p[38]))
    p32 = p[31]
    p18 = p[17]
    p14 = p[13]
    p39 = p[38]

    # F10.7
    t[0] = p[21] * dfa

    # Time Independent
    t[1] = (
        p[1] * plg[0][2]
        + p[2] * plg[0][4]
        + p[22] * plg[0][6]
        + p[26] * plg[0][1]
        + p[14] * plg[0][3]
        + p[59] * plg[0][5]
    )

    # Symmetrical Annual
    t[2] = (p[18] + p[47] * plg[0][2] + p[29] * plg[0][4]) * cd32

    # Symmetrical Semi-Annual
    t[3] = (p[15] + p[16] * plg[0][2] + p[30] * plg[0][4]) * cd18

    # Asymmetrical Annual
    t[4] = (p[9] * plg[0][1] + p[10] * plg[0][3] + p[20] * plg[0][5]) * cd14

    # Asymmetrical Semi-Annual
    t[5] = (p[37] * plg[0][1]) * cd39

    # DIURNAL
    if flags.sw[7]:
        t71 = p[11] * plg[1][2] * cd14 * flags.swc[5]
        t72 = p[12] * plg[1][2] * cd14 * flags.swc[5]
        t[6] = (p[3] * plg[1][1] + p[4] * plg[1][3] + t71) * ctloc + (
            p[6] * plg[1][1] + p[7] * plg[1][3] + t72
        ) * stloc

    # Semi-DIURNAL */
    if flags.sw[8]:
        t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * flags.swc[5]
        t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * flags.swc[5]
        t[7] = (p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc + (
            p[8] * plg[2][2] + p[42] * plg[2][4] + t82
        ) * s2tloc

    # TerDIURNAL */
    if flags.sw[14]:
        t[13] = p[39] * plg[3][3] * s3tloc + p[40] * plg[3][3] * c3tloc

    # Magnetic Activity
    if flags.sw[9]:
        if flags.sw[9] == 1:
            t[8] = apdf * (p[32] + p[45] * plg[0][2] * flags.swc[2])
        if flags.sw[9] == -1:
            t[8] = p[50] * apt[0] + p[96] * plg[0][2] * apt[0] * flags.swc[2]

    # Longitudinal
    if not ((flags.sw[10] == 0) or (flags.sw[11] == 0) or (input.g_long <= -1000.0)):
        t[10] = (
            1.0
            + plg[0][1]
            * (
                p[80] * flags.swc[5] * np.cos(dr * (input.doy - p[81]))
                + p[85] * flags.swc[6] * np.cos(2.0 * dr * (input.doy - p[86]))
            )
            + p[83] * flags.swc[3] * np.cos(dr * (input.doy - p[84]))
            + p[87] * flags.swc[4] * np.cos(2.0 * dr * (input.doy - p[88]))
        ) * (
            (
                p[64] * plg[1][2]
                + p[65] * plg[1][4]
                + p[66] * plg[1][6]
                + p[74] * plg[1][1]
                + p[75] * plg[1][3]
                + p[76] * plg[1][5]
            )
            * np.cos(dgtr * input.g_long)
            + (
                p[90] * plg[1][2]
                + p[91] * plg[1][4]
                + p[92] * plg[1][6]
                + p[77] * plg[1][1]
                + p[78] * plg[1][3]
                + p[79] * plg[1][5]
            )
            * np.sin(dgtr * input.g_long)
        )

    tt = 0
    for ind in range(0, 14):
        tt += abs(flags.sw[ind + 1]) * t[ind]
    return tt


# GTD7
def gtd7(input, flags, output):
    """
    Neutral Atmosphere Empircial Model from the surface to lower exosphere.

    Parameters
    ----------
    input : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_input
        Nrlmsise_input class instance
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Nrlmsise_flags class instance
    output :  ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_output
        Nrlmsise_output class instance
    """
    mn3 = 5
    zn3 = [32.5, 20.0, 15.0, 10.0, 0.0]
    mn2 = 4
    zn2 = [72.5, 55.0, 45.0, 32.5]
    zmix = 62.5
    soutput = nrlmsise_output()

    tselec(flags)

    # Latitude Variation of Gravity (none for sw[2] = 0)
    xlat = input.g_lat
    if flags.sw[2] == 0:
        xlat = 45.0
    glatf(xlat, gsurf, re)

    xmm = pdm[2][4]

    # Thermosphere / Mesosphere (Above ZN2[0])
    if input.alt > zn2[0]:
        altt = input.alt
    else:
        altt = zn2[0]

    tmp = input.alt
    input.alt = altt

    gts7(input, flags, soutput)

    altt = input.alt
    input.alt = tmp

    if flags.sw[0]:
        """
        Metric Adjustments
        """
        dm28m = dm28 * 1.0e6
    else:
        dm28m = dm28

    output.t[0] = soutput.t[0]
    output.t[1] = soutput.t[1]

    if input.alt >= zn2[0]:
        for ind in range(0, 9):
            output.d[ind] = soutput.d[ind]
        return

    """
    Lower Mesosphere / Upper Stratosphere (Between ZN3[0] & ZN2[0])
        Temperature at Nodes & Gradients at the End Nodes
        Inverse Temperature a Linear Function of Spherical Harmonics
    """
    meso_tgn2[0] = meso_tgn1[1]
    meso_tn2[0] = meso_tn1[4]

    meso_tn2[1] = (
        pma[0][0] * pavgm[0] / (1.0 - flags.sw[20] * glob7s(pma[0], input, flags))
    )
    meso_tn2[2] = (
        pma[1][0] * pavgm[1] / (1.0 - flags.sw[20] * glob7s(pma[1], input, flags))
    )
    meso_tn2[3] = (
        pma[2][0]
        * pavgm[2]
        / (1.0 - flags.sw[20] * flags.sw[22] * glob7s(pma[2], input, flags))
    )

    meso_tgn2[1] = (
        pavgm[8]
        * pma[9][0]
        * (1.0 + flags.sw[20] * flags.sw[22] * glob7s(pma[9], input, flags))
        * meso_tn2[3]
        * meso_tn2[3]
        / pow((pma[2][0] * pavgm[2]), 2.0)
    )
    meso_tn3[0] = meso_tn2[3]

    if input.alt < zn3[0]:
        """
        Lower Stratosphere and Troposphere (Below ZN3[0])
            Temperature at Nodes and Gradients at End Nodes
            Inverse Temperature a Linear Function of Spherical Harmonics
        """
        meso_tgn3[0] = meso_tgn2[1]
        meso_tn3[1] = (
            pma[3][0] * pavgm[3] / (1.0 - flags.sw[22] * glob7s(pma[3], input, flags))
        )
        meso_tn3[2] = (
            pma[4][0] * pavgm[4] / (1.0 - flags.sw[22] * glob7s(pma[4], input, flags))
        )
        meso_tn3[3] = (
            pma[5][0] * pavgm[5] / (1.0 - flags.sw[22] * glob7s(pma[5], input, flags))
        )
        meso_tn3[4] = (
            pma[6][0] * pavgm[6] / (1.0 - flags.sw[22] * glob7s(pma[6], input, flags))
        )
        meso_tgn3[1] = (
            pma[7][0]
            * pavgm[7]
            * (1.0 + flags.sw[22] * glob7s(pma[7], input, flags))
            * meso_tn3[4]
            * meso_tn3[4]
            / pow((pma[6][0] * pavgm[6]), 2.0)
        )

    """
    Linear Transition To Full Mixing Below ZN2[0]
    """

    dmc = 0
    if input.alt > zmix:
        dmc = 1.0 - (zn2[0] - input.alt) / (zn2[0] - zmix)
    dz28 = soutput.d[2]

    # N2 Density
    dmr = soutput.d[2] / dm28m - 1.0
    tz = [0.0]
    output.d[2] = densm(
        input.alt,
        dm28m,
        xmm,
        tz,
        mn3,
        zn3,
        meso_tn3,
        meso_tgn3,
        mn2,
        zn2,
        meso_tn2,
        meso_tgn2,
    )
    output.d[2] = output.d[2] * (1.0 + dmr * dmc)

    # HE Density
    dmr = soutput.d[0] / (dz28 * pdm[0][1]) - 1.0
    output.d[0] = output.d[2] * pdm[0][1] * (1.0 + dmr * dmc)

    # O Density
    output.d[1] = 0
    output.d[8] = 0

    # O2 Density
    dmr = soutput.d[3] / (dz28 * pdm[3][1]) - 1.0
    output.d[3] = output.d[2] * pdm[3][1] * (1.0 + dmr * dmc)

    # AR Density
    dmr = soutput.d[4] / (dz28 * pdm[4][1]) - 1.0
    output.d[4] = output.d[2] * pdm[4][1] * (1.0 + dmr * dmc)

    # HYDROGEN Density
    output.d[6] = 0

    # ATOMIC NITROGEN Density
    output.d[7] = 0

    # Total Mass Density
    output.d[5] = 1.66e-24 * (
        4.0 * output.d[0]
        + 16.0 * output.d[1]
        + 28.0 * output.d[2]
        + 32.0 * output.d[3]
        + 40.0 * output.d[4]
        + output.d[6]
        + 14.0 * output.d[7]
    )

    if flags.sw[0]:
        output.d[5] = output.d[5] / 1000

    # Temperature At Altitude
    global dd
    dd = densm(
        input.alt,
        1.0,
        0,
        tz,
        mn3,
        zn3,
        meso_tn3,
        meso_tgn3,
        mn2,
        zn2,
        meso_tn2,
        meso_tgn2,
    )
    output.t[1] = tz[0]
    return


# GTD7D
def gtd7d(input, flags, output):
    """
    This function provides Effective Total Mass Density for output d[5] 
    which includes contributions from "anomalous oxygen" which can 
    affect satellite drag above 500 km. See the section "output" for
    additional details.

    Parameters
    ----------
    input : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_input
        Nrlmsise_input class instance
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Nrlmsise_flags class instance
    output :  ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_output
        Nrlmsise_output class instance
    """
    gtd7(input, flags, output)
    output.d[5] = 1.66e-24 * (
        4.0 * output.d[0]
        + 16.0 * output.d[1]
        + 28.0 * output.d[2]
        + 32.0 * output.d[3]
        + 40.0 * output.d[4]
        + output.d[6]
        + 14.0 * output.d[7]
        + 16.0 * output.d[8]
    )
    if flags.sw[0]:
        output.d[5] = output.d[5] / 1000
    return


# GHP7
def ghp7(input, flags, output, press):
    """
    To specify outputs at a pressure level (press) rather than at 
    an altitude.

    Parameters
    ----------
    input : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_input
        Nrlmsise_input class instance
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Nrlmsise_flags class instance
    output :  ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_output
        Nrlmsise_output class instance
    press : float
    """
    bm = 1.3806e-19
    test = 0.00043
    ltest = 12

    pl = math.log10(press)

    if pl >= -5.0:
        if pl > 2.5:
            zi = 18.06 * (3.00 - pl)
        elif (pl > 0.075) and (pl <= 2.5):
            zi = 14.98 * (3.08 - pl)
        elif (pl > -1) and (pl <= 0.075):
            zi = 17.80 * (2.72 - pl)
        elif (pl > -2) and (pl <= -1):
            zi = 14.28 * (3.64 - pl)
        elif (pl > -4) and (pl <= -2):
            zi = 12.72 * (4.32 - pl)
        elif pl <= -4:
            zi = 25.3 * (0.11 - pl)

        cl = input.g_lat / 90.0
        cl2 = cl * cl

        if input.doy < 182:
            cd = (1.0 - float(input.doy)) / 91.25
        else:
            cd = (float(input.doy)) / 91.25 - 3.0

        ca = 0

        if (pl > -1.11) and (pl <= -0.23):
            ca = 1.0
        if pl > -0.23:
            ca = (2.79 - pl) / (2.79 + 0.23)
        if (pl <= -1.11) and (pl > -3):
            ca = (-2.93 - pl) / (-2.93 + 1.11)
        z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl

    else:
        z = 22.0 * pow((pl + 4.0), 2.0) + 110.0

    # Iteration Loop
    l = 0
    while True:
        l += 1
        input.alt = z

        gtd7(input, flags, output)

        z = input.alt
        xn = (
            output.d[0]
            + output.d[1]
            + output.d[2]
            + output.d[3]
            + output.d[4]
            + output.d[6]
            + output.d[7]
        )
        p = bm * xn * output.t[1]
        if flags.sw[0]:
            p = p * 1.0e-6
        diff = pl - math.log10(p)

        if math.sqrt(diff * diff) < test:
            return
        if l == ltest:
            print(f"ERROR: ghp7 not converging for press {press} diff {diff}")
            return

        xm = output.d[5] / xn / 1.66e-24
        if flags.sw[0]:
            xm = xm * 1.0e3
        g = gsurf[0] / pow((1.0 + z / re[0]), 2.0)
        sh = rgas * output.t[1] / (xm * g)

        # New Altitude Estimate Using Scale Height
        if l < 6:
            z = z - sh * diff * 2.302
        else:
            z = z - sh * diff


# GTS7
def gts7(input, flags, output):
    """
    Thermospheric Portion Of NRLMSISE-00
        Check GTD7 for more extensive comments
        alt > 72.5 KM !!

    Parameters
    ----------
    input : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_input
        Nrlmsise_input class instance
    flags : ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_flags
        Nrlmsise_flags class instance
    output :  ~poliastro.earth.atmosphere.nrlmsise00.nrlmsise_output
        Nrlmsise_output class instance
    """
    mn1 = 5
    zn1 = [120.0, 110.0, 100.0, 90.0, 72.5]
    alpha = [-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0]
    altl = [200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0]
    za = pdl[1][15]
    zn1[0] = za

    for ind in range(9):
        output.d[ind] = 0

    # TINF Variations Not Important Below ZA or ZN1(1)
    if input.alt > zn1[0]:
        tinf = ptm[0] * pt[0] * (1.0 + flags.sw[16] * globe7(pt, input, flags))
    else:
        tinf = ptm[0] * pt[0]
    output.t[0] = tinf

    # Gradient Variation Not Important Below ZN1(5)
    if input.alt > zn1[4]:
        g0 = ptm[3] * ps[0] * (1.0 + flags.sw[19] * globe7(ps, input, flags))
    else:
        g0 = ptm[3] * ps[0]
    tlb = ptm[1] * (1.0 + flags.sw[17] * globe7(pd[3], input, flags)) * pd[3][0]
    s = g0 / (tinf - tlb)

    """
    Lower Thermosphere Temperature Variations Not Significant for Density Above 300 KM
    """
    if input.alt < 300.0:
        meso_tn1[1] = (
            ptm[6] * ptl[0][0] / (1.0 - flags.sw[18] * glob7s(ptl[0], input, flags))
        )
        meso_tn1[2] = (
            ptm[2] * ptl[1][0] / (1.0 - flags.sw[18] * glob7s(ptl[1], input, flags))
        )
        meso_tn1[3] = (
            ptm[7] * ptl[2][0] / (1.0 - flags.sw[18] * glob7s(ptl[2], input, flags))
        )
        meso_tn1[4] = (
            ptm[4]
            * ptl[3][0]
            / (1.0 - flags.sw[18] * flags.sw[20] * glob7s(ptl[3], input, flags))
        )
        meso_tgn1[1] = (
            ptm[8]
            * pma[8][0]
            * (1.0 + flags.sw[18] * flags.sw[20] * glob7s(pma[8], input, flags))
            * meso_tn1[4]
            * meso_tn1[4]
            / pow((ptm[4] * ptl[3][0]), 2.0)
        )
    else:
        meso_tn1[1] = ptm[6] * ptl[0][0]
        meso_tn1[2] = ptm[2] * ptl[1][0]
        meso_tn1[3] = ptm[7] * ptl[2][0]
        meso_tn1[4] = ptm[4] * ptl[3][0]
        meso_tgn1[1] = (
            ptm[8]
            * pma[8][0]
            * meso_tn1[4]
            * meso_tn1[4]
            / pow((ptm[4] * ptl[3][0]), 2.0)
        )

    z0 = zn1[3]
    t0 = meso_tn1[3]
    tr12 = 1.0

    # N2 Variation Factor at ZLB
    g28 = flags.sw[21] * globe7(pd[2], input, flags)

    # Variation Of TURBOPAUSE Height
    zhf = pdl[1][24] * (
        1.0
        + flags.sw[5]
        * pdl[0][24]
        * np.sin(dgtr * input.g_lat)
        * np.cos(dr * (input.doy - pt[13]))
    )
    output.t[0] = tinf
    xmm = pdm[2][4]
    z = input.alt

    # N2 Density
    # Diffusive Density at ZLB
    db28 = pdm[2][0] * np.exp(g28) * pd[2][0]

    # Diffusive Density at Alt
    temp = [output.t[1]]
    output.d[2] = densu(
        z,
        db28,
        tinf,
        tlb,
        28.0,
        alpha[2],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[2]

    # Turbopause
    zh28 = pdm[2][2] * zhf
    zhm28 = pdm[2][3] * pdl[1][5]
    xmd = 28.0 - xmm

    # Mixed Density at Zlb
    tz = [0]
    b28 = densu(
        zh28,
        db28,
        tinf,
        tlb,
        xmd,
        (alpha[2] - 1.0),
        tz,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    if (flags.sw[15]) and (z <= altl[2]):
        # Mixed Density at Alt
        global dm28
        dm28 = densu(
            z,
            b28,
            tinf,
            tlb,
            xmm,
            alpha[2],
            tz,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        # Net Density at Alt
        output.d[2] = dnet(output.d[2], dm28, zhm28, xmm, 28.0)

    # HE DENSITY
    # Density Variation Factor at ZLB
    g4 = flags.sw[21] * globe7(pd[0], input, flags)

    #  Diffusive Density at ZLB
    db04 = pdm[0][0] * np.exp(g4) * pd[0][0]

    #  Diffusive Density at Alt
    temp = [output.t[1]]
    output.d[0] = densu(
        z,
        db04,
        tinf,
        tlb,
        4.0,
        alpha[0],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[0]
    if (flags.sw[15]) and (z < altl[0]):
        # Turbopause
        zh04 = pdm[0][2]
        # Mixed density at ZLB
        temp = [output.t[1]]
        b04 = densu(
            zh04,
            db04,
            tinf,
            tlb,
            4.0 - xmm,
            alpha[0] - 1.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        # Mixed density at Alt
        temp = [output.t[1]]
        global dm04
        dm04 = densu(
            z,
            b04,
            tinf,
            tlb,
            xmm,
            0.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        zhm04 = zhm28
        # Net density at Alt
        output.d[0] = dnet(output.d[0], dm04, zhm04, xmm, 4.0)
        # Correction to specified mixing ratio at ground
        rl = math.log(b28 * pdm[0][1] / b04)
        zc04 = pdm[0][4] * pdl[1][0]
        hc04 = pdm[0][5] * pdl[1][1]
        # Net density corrected at Alt
        output.d[0] = output.d[0] * ccor(z, rl, hc04, zc04)

    # # O DENSITY
    # # Density Variation Factor at ZLB
    g16 = flags.sw[21] * globe7(pd[1], input, flags)

    # Diffusive Density at ZLB
    db16 = pdm[1][0] * np.exp(g16) * pd[1][0]

    # Diffusive density at Alt
    temp = [output.t[1]]
    output.d[1] = densu(
        z,
        db16,
        tinf,
        tlb,
        16.0,
        alpha[1],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[1]
    if (flags.sw[15]) and (z <= altl[1]):
        # Turbopause
        zh16 = pdm[1][2]
        # Mixed density at ZLB
        temp = [output.t[1]]
        b16 = densu(
            zh16,
            db16,
            tinf,
            tlb,
            16.0 - xmm,
            (alpha[1] - 1.0),
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        # Mixed density at Alt
        temp = [output.t[1]]
        global dm16
        dm16 = densu(
            z,
            b16,
            tinf,
            tlb,
            xmm,
            0.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        zhm16 = zhm28
        # Net density at Alt
        output.d[1] = dnet(output.d[1], dm16, zhm16, xmm, 16.0)
        rl = (
            pdm[1][1]
            * pdl[1][16]
            * (1.0 + flags.sw[1] * pdl[0][23] * (input.f107A - 150.0))
        )
        hc16 = pdm[1][5] * pdl[1][3]
        zc16 = pdm[1][4] * pdl[1][2]
        hc216 = pdm[1][5] * pdl[1][4]
        output.d[1] = output.d[1] * ccor2(z, rl, hc16, zc16, hc216)
        # Chemistry correction
        hcc16 = pdm[1][7] * pdl[1][13]
        zcc16 = pdm[1][6] * pdl[1][12]
        rc16 = pdm[1][3] * pdl[1][14]
        # Net density corrected at Alt
        output.d[1] = output.d[1] * ccor(z, rc16, hcc16, zcc16)

    # O2 DENSITY
    # Density Variation Factor at ZLB
    g32 = flags.sw[21] * globe7(pd[4], input, flags)

    # Diffusive Density at ZLB
    db32 = pdm[3][0] * np.exp(g32) * pd[4][0]

    # Diffusive Density at ALT
    temp = [output.t[1]]
    output.d[3] = densu(
        z,
        db32,
        tinf,
        tlb,
        32.0,
        alpha[3],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[3]
    if flags.sw[15]:
        if z <= altl[3]:
            # Turbopause
            zh32 = pdm[3][2]
            # Mixed Density at ZLB
            temp = [output.t[1]]
            b32 = densu(
                zh32,
                db32,
                tinf,
                tlb,
                32.0 - xmm,
                alpha[3] - 1.0,
                temp,
                ptm[5],
                s,
                mn1,
                zn1,
                meso_tn1,
                meso_tgn1,
            )
            output.t[1] = temp[0]
            # Mixed Density at ALT
            temp = [output.t[1]]
            global dm32
            dm32 = densu(
                z,
                b32,
                tinf,
                tlb,
                xmm,
                0.0,
                temp,
                ptm[5],
                s,
                mn1,
                zn1,
                meso_tn1,
                meso_tgn1,
            )
            output.t[1] = temp[0]
            zhm32 = zhm28
            # Net Density at ALT
            output.d[3] = dnet(output.d[3], dm32, zhm32, xmm, 32.0)
            # Correction to Specific Mixing Ratio at Ground
            rl = math.log(b28 * pdm[3][1] / b32)
            hc32 = pdm[3][5] * pdl[1][7]
            zc32 = pdm[3][4] * pdl[1][6]
            output.d[3] = output.d[3] * ccor(z, rl, hc32, zc32)

        # Correction for General Departure From Diffusive Equilibrium Above ZLB
        hcc32 = pdm[3][7] * pdl[1][22]
        hcc232 = pdm[3][7] * pdl[0][22]
        zcc32 = pdm[3][6] * pdl[1][21]
        rc32 = (
            pdm[3][3]
            * pdl[1][23]
            * (1.0 + flags.sw[1] * pdl[0][23] * (input.f107A - 150.0))
        )

        # Net Density Corrected at ALT
        output.d[3] = output.d[3] * ccor2(z, rc32, hcc32, zcc32, hcc232)

    # AR DENSITY
    # Density Variation Factor at ZLB
    g40 = flags.sw[21] * globe7(pd[5], input, flags)

    # Diffusive Density at ZLB
    db40 = pdm[4][0] * np.exp(g40) * pd[5][0]

    # Diffusive Density at ALT
    temp = [output.t[1]]
    output.d[4] = densu(
        z,
        db40,
        tinf,
        tlb,
        40.0,
        alpha[4],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[4]
    if (flags.sw[15]) and (z <= altl[4]):
        # Turbopause
        zh40 = pdm[4][2]
        # Mixed Density at ZLB
        temp = [output.t[1]]
        b40 = densu(
            zh40,
            db40,
            tinf,
            tlb,
            40.0 - xmm,
            alpha[4] - 1.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        # Mixed Density at ALT
        temp = [output.t[1]]
        global dm40
        dm40 = densu(
            z,
            b40,
            tinf,
            tlb,
            xmm,
            0.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        zhm40 = zhm28
        # Net Density at ALT
        output.d[4] = dnet(output.d[4], dm40, zhm40, xmm, 40.0)
        # Correction To Specified Mixing Ratio at Ground
        rl = math.log(b28 * pdm[4][1] / b40)
        hc40 = pdm[4][5] * pdl[1][9]
        zc40 = pdm[4][4] * pdl[1][8]
        # Net Density Corrected at ALT
        output.d[4] = output.d[4] * ccor(z, rl, hc40, zc40)

    # HYDROGEN DENSITY
    # Density Variation Factor at ZLB
    g1 = flags.sw[21] * globe7(pd[6], input, flags)

    # Diffusive Density at ZLB
    db01 = pdm[5][0] * np.exp(g1) * pd[6][0]

    # Diffusive Density at ALT
    temp = [output.t[1]]
    output.d[6] = densu(
        z,
        db01,
        tinf,
        tlb,
        1.0,
        alpha[6],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[6]
    if (flags.sw[15]) and (z <= altl[6]):
        # Turbopause
        zh01 = pdm[5][2]
        # Mixed Density at ZLB
        temp = [output.t[1]]
        b01 = densu(
            zh01,
            db01,
            tinf,
            tlb,
            1.0 - xmm,
            alpha[6] - 1.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        # Mixed Density at ALT
        temp = [output.t[1]]
        global dm01
        dm01 = densu(
            z,
            b01,
            tinf,
            tlb,
            xmm,
            0.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        zhm01 = zhm28
        # Net Density at ALT
        output.d[6] = dnet(output.d[6], dm01, zhm01, xmm, 1.0)
        # Correction to Specified Mixing Ratio at
        rl = math.log(b28 * pdm[5][1] * math.sqrt(pdl[1][17] * pdl[1][17]) / b01)
        hc01 = pdm[5][5] * pdl[1][11]
        zc01 = pdm[5][4] * pdl[1][10]
        output.d[6] = output.d[6] * ccor(z, rl, hc01, zc01)
        # Chemistry Correction
        hcc01 = pdm[5][7] * pdl[1][19]
        zcc01 = pdm[5][6] * pdl[1][18]
        rc01 = pdm[5][3] * pdl[1][20]
        # Net Density Corrected at ALT
        output.d[6] = output.d[6] * ccor(z, rc01, hcc01, zcc01)

    # ATOMIC NITROGEN DENSITY
    # Density Variation Factor at ZLB
    g14 = flags.sw[21] * globe7(pd[7], input, flags)

    # Diffusive Density at ZLB
    db14 = pdm[6][0] * np.exp(g14) * pd[7][0]

    # Diffusive Density at ALT
    temp = [output.t[1]]
    output.d[7] = densu(
        z,
        db14,
        tinf,
        tlb,
        14.0,
        alpha[7],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    dd = output.d[7]
    if (flags.sw[15]) and (z <= altl[7]):
        # Turbopause
        zh14 = pdm[6][2]
        # Mixed Density at ZLB
        temp = [output.t[1]]
        b14 = densu(
            zh14,
            db14,
            tinf,
            tlb,
            14.0 - xmm,
            alpha[7] - 1.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        # Mixed Density at ALT
        temp = [output.t[1]]
        global dm14
        dm14 = densu(
            z,
            b14,
            tinf,
            tlb,
            xmm,
            0.0,
            temp,
            ptm[5],
            s,
            mn1,
            zn1,
            meso_tn1,
            meso_tgn1,
        )
        output.t[1] = temp[0]
        zhm14 = zhm28
        # Net Density at ALT
        output.d[7] = dnet(output.d[7], dm14, zhm14, xmm, 14.0)
        # Correction to Specified Mixing Ratio at Ground
        rl = math.log(b28 * pdm[6][1] * math.sqrt(pdl[0][2] * pdl[0][2]) / b14)
        hc14 = pdm[6][5] * pdl[0][1]
        zc14 = pdm[6][4] * pdl[0][0]
        output.d[7] = output.d[7] * ccor(z, rl, hc14, zc14)
        # Chemistry Correction
        hcc14 = pdm[6][7] * pdl[0][4]
        zcc14 = pdm[6][6] * pdl[0][3]
        rc14 = pdm[6][3] * pdl[0][5]
        # Net Density Corrected at ALT
        output.d[7] = output.d[7] * ccor(z, rc14, hcc14, zcc14)

    # ANOMALOUS OXYGEN DENSITY
    g16h = flags.sw[21] * globe7(pd[8], input, flags)
    db16h = pdm[7][0] * np.exp(g16h) * pd[8][0]
    tho = pdm[7][9] * pdl[0][6]
    temp = [output.t[1]]
    dd = densu(
        z,
        db16h,
        tho,
        tho,
        16.0,
        alpha[8],
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    zsht = pdm[7][5]
    zmho = pdm[7][4]
    zsho = scalh(zmho, 16.0, tho)
    output.d[8] = dd * np.exp(-zsht / zsho * (np.exp(-(z - zmho) / zsht) - 1.0))

    # Total Mass Density
    output.d[5] = 1.66e-24 * (
        4.0 * output.d[0]
        + 16.0 * output.d[1]
        + 28.0 * output.d[2]
        + 32.0 * output.d[3]
        + 40.0 * output.d[4]
        + output.d[6]
        + 14.0 * output.d[7]
    )
    db48 = 1.66e-24 * (
        4.0 * db04
        + 16.0 * db16
        + 28.0 * db28
        + 32.0 * db32
        + 40.0 * db40
        + db01
        + 14.0 * db14
    )

    # Temperature
    z = math.sqrt(input.alt * input.alt)
    temp = [output.t[1]]
    ddum = densu(
        z,
        1.0,
        tinf,
        tlb,
        0.0,
        0.0,
        temp,
        ptm[5],
        s,
        mn1,
        zn1,
        meso_tn1,
        meso_tgn1,
    )
    output.t[1] = temp[0]
    if flags.sw[0]:
        for ind in range(0, 9):
            output.d[ind] = output.d[ind] * 1.0e6
        output.d[5] = output.d[5] / 1000
    return
