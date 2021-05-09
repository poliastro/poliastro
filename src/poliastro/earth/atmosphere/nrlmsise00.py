## Implementation for NRLMSISE 2001 Atmosphere Model

import math

import astropy
import numpy as np
from astropy import units as u

# Constant Values (Units to be added later)
dgtr = 1.74533e-2
rgas = 831.4
sr = 7.2722e-5
dr = 1.72142e-2
hr = 0.2618


class nrlmsise_flags:
    """
    Switches: to turn on and off particular variations use these switches.
    0 is off, 1 is on, and 2 is main effects off but cross terms on.

    Standard values are 0 for switch 0 and 1 for switches 1 to 23. The
    array "switches" needs to be set accordingly by the calling program.
    The arrays sw and swc are set internally.

    switches[i]:
    i  - explanation
    -----------------
    0  - output in meters and kilograms instead of centimeters and grams
    1  - F10.7 effect on mean
    2  - time independent
    3  - symmetrical annual
    4  - symmetrical semiannual
    5  - asymmetrical annual
    6  - asymmetrical semiannual
    7  - diurnal
    8  - semidiurnal
    9  - daily ap [when this is set to -1 (!) the pointer ap_a in struct nrlmsise_input must
                    point to a struct ap_array]
    10 - all UT/long effects
    11 - longitudinal
    12 - UT and mixed UT/long
    13 - mixed AP/UT/LONG
    14 - terdiurnal
    15 - departures from diffusive equilibrium
    16 - all TINF var
    17 - all TLB var
    18 - all TN1 var
    19 - all S var
    20 - all TN2 var
    21 - all NLB var
    22 - all TN3 var
    23 - turbo scale height var
    """

    def __init__(self):
        self.switches = np.empty(24, dtype=np.int16)
        self.sw = np.empty(24, dtype=np.float16)
        self.swc = np.empty(24, dtype=np.float16)


class ap_array:
    """
    Array containing the following magnetic values:
    0 : daily AP
    1 : 3 hr AP index for current time
    2 : 3 hr AP index for 3 hrs before current time
    3 : 3 hr AP index for 6 hrs before current time
    4 : 3 hr AP index for 9 hrs before current time
    5 : Average of eight 3 hr AP indicies from 12 to 33 hrs
        prior to current time
    6 : Average of eight 3 hr AP indicies from 36 to 57 hrs
        prior to current time
    """

    def __init__(self):
        a = np.empty(7, dtype=np.float16)


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
    """

    def __init__(
        self,
        year,
        doy,
        sec,
        alt,
        g_lat,
        g_long,
        lst,
        f107A=150,
        f107=150,
        ap=4,
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
        self.ap_a = ap_array()


class nrlmsise_output:
    """
    OUTPUT VARIABLES:
    d[0] - HE NUMBER DENSITY(CM-3)
    d[1] - O NUMBER DENSITY(CM-3)
    d[2] - N2 NUMBER DENSITY(CM-3)
    d[3] - O2 NUMBER DENSITY(CM-3)
    d[4] - AR NUMBER DENSITY(CM-3)
    d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
    d[6] - H NUMBER DENSITY(CM-3)
    d[7] - N NUMBER DENSITY(CM-3)
    d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)

    t[0] - EXOSPHERIC TEMPERATURE
    t[1] - TEMPERATURE AT ALT

    O, H, and N are set to zero below 72.5 km

    t[0], Exospheric temperature, is set to global average for
    altitudes below 120 km. The 120 km gradient is left at global
    average value for altitudes below 72 km.

    d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
    and GTD7D

    SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
    species labeled by indices 0-4 and 6-7 in output variable d.
    This includes He, O, N2, O2, Ar, H, and N but does NOT include
    anomalous oxygen (species index 8).

    SUBROUTINE GTD7D -- d[5] is the "effective total mass density
    for drag" and is the sum of the mass densities of all species
    in this model, INCLUDING anomalous oxygen.
    """

    def __init__(self):
        d = np.empty(9, dtype=np.float16)
        t = np.empty(2, dtype=np.float16)


# # GTD7
# # Neutral Atmosphere Empircial Model from the surface to lower exosphere.
# def gtd7():
#     nrlmsise_input():

# void gtd7 (struct nrlmsise_input *input, \
#            struct nrlmsise_flags *flags, \
#            struct nrlmsise_output *output);


# /* GTD7D */
# /*   This subroutine provides Effective Total Mass Density for output
#  *   d[5] which includes contributions from "anomalous oxygen" which can
#  *   affect satellite drag above 500 km. See the section "output" for
#  *   additional details.
#  */
# void gtd7d(struct nrlmsise_input *input, \
#            struct nrlmsise_flags *flags, \
#            struct nrlmsise_output *output);


# /* GTS7 */
# /*   Thermospheric portion of NRLMSISE-00
#  */
# void gts7 (struct nrlmsise_input *input, \
# 	   struct nrlmsise_flags *flags, \
# 	   struct nrlmsise_output *output);


# /* GHP7 */
# /*   To specify outputs at a pressure level (press) rather than at
#  *   an altitude.
#  */
# void ghp7 (struct nrlmsise_input *input, \
#            struct nrlmsise_flags *flags, \
#            struct nrlmsise_output *output, \
#            double press);
###################################################################################
flags = nrlmsise_flags()
input = nrlmsise_input()
output = nrlmsise_output()

# MESO7
meso_tn1 = np.empty(5, dtype=np.float16)
meso_tn2 = np.empty(4, dtype=np.float16)
meso_tn3 = np.empty(5, dtype=np.float16)
meso_tgn1 = np.empty(2, dtype=np.float16)
meso_tgn2 = np.empty(2, dtype=np.float16)
meso_tgn3 = np.empty(2, dtype=np.float16)

# POWER7
pt = np.empty(150, dtype=np.float16)
pd = np.empty([9, 150], dtype=np.float16)
ps = np.empty(150, dtype=np.float16)
pdl = np.empty([2, 25], dtype=np.float16)
ptl = np.empty([4, 100], dtype=np.float16)
pma = np.empty([10, 100], dtype=np.float16)
sam = np.empty(100, dtype=np.float16)

# LOWER7
ptm = np.empty(10, dtype=np.float16)
pdm = np.empty([8, 10], dtype=np.float16)
pavgm = np.empty(10, dtype=np.float16)

# LPOLY
dfa = 0
plg = np.empty([4, 9], dtype=np.float16)
ctloc, stloc = 0
c2tloc, s2tloc = 0
s3tloc, c3tloc = 0
apdf = 0
apt = np.empty(4, dtype=np.float16)


def tselec(flags):
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


def glatf(lat, gv, reff):
    c2 = np.cos(2.0 * dgtr * lat)
    gv = 980.616 * (1.0 - 0.0026373 * c2)
    reff = 2.0 * gv / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5


def ccor(alt, r, h1, zh):
    """
    Chemistry/Dissociation Correction for MSIS Models
    alt - altitude
    r   - target ratio
    h1  - transition scale length
    zh  - altitude of (1/2 r)
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
    alt - altitude
    r   - target ratio
    h1  - transition scale length
    zh  - altitude of (1/2 r)
    h2  - transition scale length #2 ?
    """
    e1 = (alt - zh) / h1
    e2 = (alt - zh) / h2
    if e1 > 70 or e2 > 70:
        return np.exp(0)
    if e1 < -70 or e2 < -70:
        return np.exp(r)
    ex1 = np.exp(e1)
    ex2 = np.exp(e2)
    ccor2v = r / (1.0 + 0.5 * (ex1 + ex2))
    return np.exp(ccor2v)


def scalh(alt, xm, temp):
    # gsurf & re to be declared
    g = gsurf / (1.0 + (alt / re)) ** 2
    g = rgas * temp / (g * xm)
    return g


def dnet(dd, dm, zhm, xmm, xm):
    """
    TurboPause Correction For MSIS Model
    Root Mean Density
    dd   - diffusive density
    dm   - full mixed density
    zhm  - transition scale length
    xmm  - full mixed molecular weight
    xm   - species molecular weight
    dnet - combined density
    """
    a = zhm / (xmm - xm)
    if not (dm > 0 and dd > 0):
        # print(f"dnet log error {dm}, {dd}, {xm} \n")
        if dd == 0 and dm == 0:
            dd = 1
        if dm == 0:
            return dd
        if dd == 0:
            return dm
    ylog = a * math.log(dm, dd)
    if ylog < -10:
        return dd
    if ylog > 10:
        return dm
    a = dd * ((1.0 + np.exp(ylog)) ** (1.0 / a))
    return a


def splini(xa, ya, y2a, n, x, y):
    """
    Integrate Cubic Spline Function from XA(1) to X
    xa, ya  - array of tabulated functions in ascending order by x
    y2a     - array of second derivatives
    n       - size of arrays xa, ya, y2a
    x       - abscissa endpoint for integration
    y       - output value
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
        yi += (
            (1.0 - a ** 2) * ya[klo] / 2.0
            + b ** 2 * ya[khi] / 2.0
            + (
                (-(1.0 + a ** 4) / 4.0 + a ** 2 / 2.0) * y2a[klo]
                + (b ** 4 / 4.0 - b ** 2 / 2.0) * y2a[khi]
            )
            * h
            * h
            / 6.0
        ) * h
        klo += 1
        khi += 1
    y = yi


def splint(xa, ya, y2a, n, x, y):
    """
    Calculate Cubic Spline Interp Value
    Adapted From Numerical Recipes by PRESS ET AL.
    xa, ya  - array of tabulated functions in ascending order by x
    y2a     - array of second derivatives
    n       - size of arrays xa, ya, y2a
    x       - abscissa for interpolation
    y       - output value
    """
    klo = 0
    khi = n - 1
    while (khi - klo) > 1:
        k = (khi + klo) / 2
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
        + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * h * h / 6.0
    )
    y = yi


def spline(x, y, n, yp1, ypn, y2):
    """
    Calculate 2nd Derivatives of Cubic Spline Interp Function
    Adapted from Numerical Recipes by PRESS ET AL.
    x,y         - arrays of tabulated function in Ascending Order by x
    n           - sizes of arrays x,y
    yp1, ypn    - specified derivatives at X[0] and X[n-1],
                  Values >= 1E30 Signal Signal Second Derivative Zero
    y2          - output array of second derivatives
    """
    u = np.empty(n, dtype=np.float16)
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
            / (x[index] - x[index + 1])
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
    for index in range(0, n - 3, -1):
        y2[index] = y2[index] * y2[index + 1] + u[index]


def zeta(zz, zl):
    return (zz - zl) * (re + zl) / (re + zz)


def densm(alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2):
    """
    Calculate Temperature & Density Profiles for Lower Atoms.
    """
    xs = np.empty(10, dtype=np.float16)
    ys = np.empty(10, dtype=np.float16)
    y2out = np.empty(10, dtype=np.float16)
    y, yi, densm_temp = 0
    densm = d0
    if alt > zn2[0]:
        if xm == 0.0:
            return tz
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
    yd1 = -tgn2[0] / (t1 * t2) * zgdif
    yd2 = -tgn2[1] / (t2 * t2) * zgdif * ((re + z2) / (re + z1) ** 2)

    # Calculate Spline Coefficients
    spline(xs, ys, mn, yd1, yd2, y2out)
    x = zg / zgdif
    splint(xs, ys, y2out, mn, x, y)

    # Temperature at Altitude
    tz = 1.0 / y
    if not xm == 0.0:
        # Calculates Stratosphere / Mesosphere Density
        glb = gsurf / ((1.0 + z1 / re) ** 2)
        gamm = xm * glb * zgdif / rgas

        # Integrate Temperate Profile
        splini(xs, ys, y2out, mn, x, yi)
        expl = gamm * yi
        if expl > 50.0:
            expl = 50.0

        # Density at Altitude
        densm_temp = densm_temp * (t1 / tz) * np.exp(-expl)

    if alt > zn3[0]:
        if xm == 0.0:
            return tz
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

    for i in range(0, mn):
        xs[i] = zeta(zn3[i], z1) / zgdif
        ys[i] = 1.0 / tn3[i]
    yd1 = -tgn3[0] / (t1 * t2) * zgdif
    yd2 = -tgn3[1] / (t2 * t2) * zgdif * ((re + z2) / (re + z1)) ** 2

    # Calculate Spline Coefficients
    spline(xs, ys, mn, yd1, yd2, y2out)
    x = zg / zgdif
    splint(xs, ys, y2out, mn, x, y)

    # Temperature at Altitude
    tz = 1.0 / y
    if not xm == 0.0:
        # Calculate Tropospheric / Stratospheric Density
        glb = gsurf / ((1.0 + z1 / re) ** 2)
        gamm = xm * glb * zgdif / rgas

        # Integrate Temperature Profile
        splini(xs, ys, y2out, mn, x, yi)
        expl = gamm * yi
        if expl > 50.0:
            expl = 50.0

        # Density at Altitude
        densm_temp = densm_temp * (t1 / tz) * np.exp(-expl)

    if xm == 0.0:
        return tz
    else:
        return densm_temp


def densu(alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, mn1, zn1, tn1, tgn1):
    """
    Calculate Temperature & Density Profile for MSIS Model
    """
    # New Lower Thermo Polynomial
    x, y, yi = 0
    densu_temp = 1.0
    z1, t1, zgdif = 0
    mn = 0
    xs = np.empty(5, dtype=np.float16)
    ys = np.empty(5, dtype=np.float16)
    y2out = np.empty(5, dtype=np.float16)

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
    tz = tt
    densu_temp = tz

    if alt < za:
        # Calculates Temperature Below ZA
        # Temperature Gradient at ZA from Bates Profile.
        dta = (tinf - ta) * s2 * ((re + zlb) / (re + za)) ** 2
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
        yd2 = -tgn1[1] / (t2 * t2) * zgdif((re + z2) / (re + z1)) ** 2

        # Calculates Spline Coefficients
        spline(xs, ys, mn, yd1, yd2, y2out)
        x = zg / zgdif
        splint(xs, ys, y2out, mn, x, y)

        # Temperature at Altitude
        tz = 1.0 / y
        densu_temp = tz

    if xm == 0:
        return densu_temp

    # Calculate Density Above ZA
    glb = gsurf / ((1.0 + zlb / re) ** 2)
    gamma = xm * glb / (s2 * rgas * tinf)
    expl = np.exp(-s2 * gamma * zg2)
    if expl > 50.0:
        expl = 50.0
    if tt <= 0:
        expl = 50.0

    # Denstiy at Altitude
    densa = dlb * ((tlb / tt) ** (1.0 + alpha + gamma)) * expl
    densu_temp = densa
    if alt >= za:
        return densu_temp

    # Calculate Density Below ZA
    glb = gsurf / (1.0 + z1 / re) ** 2.0
    gamm = xm * glb * zgdif / rgas

    # Integrate Spline Temperatures
    splini(xs, ys, y2out, mn, x, yi)
    expl = gamm * yi
    if expl > 50.0:
        expl = 50.0
    if tz <= 0:
        expl = 50.0

    # Density at Altitude
    densu_temp = densu_temp * ((t1 / tz) ** (1.0 + alpha)) * np.exp(-expl)
    return densu_temp


# GLOBE7
# 3 hr Magnetic Activity Function
# Eq. A24d
def g0(a, p):
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
    return 1.0 + (1.0 - (ex ** 19.0)) / (1.0 - ex) * (ex ** 0.5)


# Eq. A24a
def sg0(ex, p, ap):
    return (
        g0(ap[1], p)
        + (
            g0(ap[2], p) * ex
            + g0(ap[3], p) * ex * ex
            + g0(ap[4], p) * (ex ** 3)
            + (g0(ap[5], p) * (ex ** 4) + g0(ap[6], p) * (ex ** 12))
            * (1.0 - (ex ** 8))
            / (1.0 - ex)
        )
    ) / sumex(ex)


def globe7(p, input, flags):
    """
    Calculate G(L) Function
    """
    # Upper Thermosphere Parameters
    t = np.empty(15, dtype=np.float16)
    ap = ap_array()
    tloc = input.lst

    for index in range(0, 14):
        t[index] = 0

    # Calculate Legendre Polynomials
    c = np.sin(input.g_lat * dgtr)
    s = np.cos(input.g_lat * dgtr)
    c2 = c ** 2
    c4 = c ** 4
    s2 = s ** 2

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
        stloc = np.sin(hr * tloc)
        ctloc = np.cos(hr * tloc)
        s2tloc = np.sin(2.0 * hr * tloc)
        c2tloc = np.cos(2.0 * hr * tloc)
        s3tloc = np.sin(3.0 * hr * tloc)
        c3tloc = np.cos(3.0 * hr * tloc)

    cd32 = np.cos(dr * input.doy - p[31])
    cd18 = np.cos(2.0 * dr * (input.doy - p[17]))
    cd14 = np.cos(dr * input.doy - p[13])
    cd39 = np.cos(2.0 * dr * (input.doy - p[38]))

    # F10.7 EFFECT
    df = input.f107 - input.f107A
    dfa = input.f107A - 150.0
    t[0] = (
        p[19] * df * (1.0 + p[59] * dfa)
        + p[20] * df * df
        + p[21] * dfa
        + p[29] * (dfa ** 2)
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
    t[2] = (p[15] + p[16] * plg[0][2]) * cd18

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
        if not p[51] == 0:
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
            p44=p[43]
            p45=p[44]
            if p44<0: