Dastcom5 parameters
===================

**avail**:

* *a* if it is available for asteroids.
* *c* if it is available for comets.
* *[number]+* since which version of DASTCOM5 is available.


+------+--------+-----------------------------------------------------------+
| avail|  Label | Definition                                                |
+======+========+===========================================================+
| ac/3+|   EPOCH|  Time of osc. orbital elements solution, JD (CT,TDB)      |
+------+--------+-----------------------------------------------------------+
| ac/3+|  CALEPO|  Time of osc. orbital elements solution, YYYYDDMM.ffff    |
+------+--------+-----------------------------------------------------------+
| ac/3+|      MA|  Mean anomaly at EPOCH, deg (elliptical & hyperbolic cases|
|      |        |  "9.999999E99" if not available)                          |
+------+--------+-----------------------------------------------------------+
| ac/3+|       W|  Argument of periapsis at EPOCH, J2000 ecliptic, deg.     |
+------+--------+-----------------------------------------------------------+
| ac/3+|      OM|  Longitude of ascending node at EPOCH, J2000 ecliptic,deg.|
+------+--------+-----------------------------------------------------------+
| ac/3+|      IN|  Inclination angle at EPOCH wrt J2000 ecliptic, deg.      |
+------+--------+-----------------------------------------------------------+
| ac/3+|      EC|  Eccentricity at EPOCH                                    |
+------+--------+-----------------------------------------------------------+
| ac/3+|       A|  Semi-major axis at EPOCH, au                             |
+------+--------+-----------------------------------------------------------+
| ac/3+|      QR|  Perihelion distance at EPOCH, au                         |
+------+--------+-----------------------------------------------------------+
| ac/3+|      TP|  Perihelion date for QR at EPOCH, JD (CT,TDB)             |
+------+--------+-----------------------------------------------------------+
| ac/3+|   TPCAL|  Perihelion date for QR at EPOCH, format YYYYMMDD.fff     |
+------+--------+-----------------------------------------------------------+
| ac/5+|  TPFRAC|  Decimal (fractional) part of TP for extended precision   |
+------+--------+-----------------------------------------------------------+
| ac/4+|  SOLDAT|  Date orbit solution was performed, JD (CT,TDB)           |
+------+--------+-----------------------------------------------------------+
| ac/4+| SRC(01)|  Square root covariance vector. Vector-stored upper-      |
|      |        |  triangular matrix with order {EC,QR,TP,OM,W,IN,{ ESTL }} |
+------+--------+-----------------------------------------------------------+
| ac/3+|       H|  Absolute visual magnitude (IAU H-G system) (99=unknown)  |
+------+--------+-----------------------------------------------------------+
| ac/3+|       G|  Mag. slope parm. (IAU H-G)(99=unknown & 0.15 not assumed)|
+------+--------+-----------------------------------------------------------+
| c/3+ |      M1|  Total absolute magnitude, mag.                           |
+------+--------+-----------------------------------------------------------+
| c/3+ |      M2|  Nuclear absolute magnitue, mag.                          |
+------+--------+-----------------------------------------------------------+
| c/4+ |      K1|  Total absolute magnitude scaling factor                  |
+------+--------+-----------------------------------------------------------+
| c/4+ |      K2|  Nuclear absolute magnitude scaling factor                |
+------+--------+-----------------------------------------------------------+
| c/4+ |   PHCOF|  Phase coefficient for K2= 5                              |
+------+--------+-----------------------------------------------------------+
| ac/3+|      A1|  Non-grav. accel., radial component, [s:10^-8 au/day^2]   |
+------+--------+-----------------------------------------------------------+
| ac/3+|      A2|  Non-grav. accel., transverse component,[s:10^-8 au/day^2 |
+------+--------+-----------------------------------------------------------+
| ac/4+|      A3|  Non-grav. accel., normal component, [s:10^-8 au/day^2]   |
+------+--------+-----------------------------------------------------------+
| c/4+ |      DT|  Non-grav. lag/delay parameter, days                      |
+------+--------+-----------------------------------------------------------+
| ac/5+|      R0|  Non-grav. model constant, normalizing distance, au       |
+------+--------+-----------------------------------------------------------+
| ac/5+|     ALN|  Non-grav. model constant, normalizing factor             |
+------+--------+-----------------------------------------------------------+
| ac/5+|      NM|  Non-grav. model constant, exponent m                     |
+------+--------+-----------------------------------------------------------+
| ac/5+|      NN|  Non-grav. model constant, exponent n                     |
+------+--------+-----------------------------------------------------------+
| ac/5+|      NK|  Non-grav. model constant, exponent k                     |
+------+--------+-----------------------------------------------------------+
| c/4+ |      S0|  Center-of-light estimated offset at 1 au, km             |
+------+--------+-----------------------------------------------------------+
| c/5+ |     TCL|  Center-of-light start-time offset, d since "ref.time"    |
+------+--------+-----------------------------------------------------------+
| a /5+|     LGK|  Surface thermal conductivity log_10(k), (W/m/K)          |
+------+--------+-----------------------------------------------------------+
| ac/5+|     RHO|  Bulk density, kg/m^3                                     |
+------+--------+-----------------------------------------------------------+
| ac/5+|   AMRAT|  Solar pressure model, area/mass ratio, m^2/kg            |
+------+--------+-----------------------------------------------------------+
| c/5+ |     AJ1|  Jet 1 acceleration, au/d^2                               |
+------+--------+-----------------------------------------------------------+
| c/5+ |     AJ2|  Jet 2 acceleration, au/d^2                               |
+------+--------+-----------------------------------------------------------+
| c/5+ |     ET1|  Thrust angle, colatitude of jet 1, deg.                  |
+------+--------+-----------------------------------------------------------+
| c/5+ |     ET2|  Thrust angle, colatitude of jet 2, deg.                  |
+------+--------+-----------------------------------------------------------+
| c/5+ |     DTH|  Jet model diurnal lag angle, deg. (delta_theta)          |
+------+--------+-----------------------------------------------------------+
| ac/5+|     ALF|  Spin pole orientation, RA, deg.                          |
+------+--------+-----------------------------------------------------------+
| ac/5+|     DEL|  Spin pole orientation, DEC, deg.                         |
+------+--------+-----------------------------------------------------------+
| ac/5+|  SPHLM3|  Earth gravity sph. harm. model limit, Earth radii        |
+------+--------+-----------------------------------------------------------+
| ac/5+|  SPHLM5|  Jupiter grav. sph. harm. model limit, Jupiter radii      |
+------+--------+-----------------------------------------------------------+
| ac/3+|      RP|  Object rotational period, hrs                            |
+------+--------+-----------------------------------------------------------+
| ac/3+|      GM|  Object mass parameter, km^3/s^2                          |
+------+--------+-----------------------------------------------------------+
| ac/3+|     RAD|  Object mean radius, km                                   |
+------+--------+-----------------------------------------------------------+
| ac/5+|  EXTNT1|  Triaxial ellipsoid, axis 1/largest equat. extent, km     |
+------+--------+-----------------------------------------------------------+
| ac/5+|  EXTNT2|  Triaxial ellipsoid, axis 2/smallest equat. extent, km    |
+------+--------+-----------------------------------------------------------+
| ac/5+|  EXTNT3|  Triaxial ellipsoid, axis 3/polar extent, km              |
+------+--------+-----------------------------------------------------------+
| ac/4+|    MOID|  Earth MOID at EPOCH time, au; '99' if not computed       |
+------+--------+-----------------------------------------------------------+
| ac/3+|  ALBEDO|  Geometric visual albedo, 99 if unknown                   |
+------+--------+-----------------------------------------------------------+
| a /3+|    BVCI|  B-V color index, mag., 99 if unknown                     |
+------+--------+-----------------------------------------------------------+
| a /5+|    UBCI|  U-B color index, mag., 99 if unknown                     |
+------+--------+-----------------------------------------------------------+
| a /5+|    IRCI|  I-R color index, mag., 99 if unknown                     |
+------+--------+-----------------------------------------------------------+
| ac/4+|    RMSW|  RMS of weighted optical residuals, arcsec                |
+------+--------+-----------------------------------------------------------+
| ac/5+|    RMSU|  RMS of unweighted optical residuals, arcsec              |
+------+--------+-----------------------------------------------------------+
| ac/5+|    RMSN|  RMS of normalized optical residuals                      |
+------+--------+-----------------------------------------------------------+
| ac/5+|   RMSNT|  RMS of all normalized residuals                          |
+------+--------+-----------------------------------------------------------+
| a /5+|    RMSH|  RMS of abs. visual magnitude (H) residuals, mag.         |
+------+--------+-----------------------------------------------------------+
| c/5+ |   RMSMT|  RMS of MT estimate residuals, mag.                       |
+------+--------+-----------------------------------------------------------+
| c/5+ |   RMSMN|  RMS of MN estimate residuals, mag.                       |
+------+--------+-----------------------------------------------------------+
| ac/3+|      NO|  Logical record-number of this object in DASTCOM          |
+------+--------+-----------------------------------------------------------+
| ac/4+|    NOBS|  Number of observations of all types used in orbit soln.  |
+------+--------+-----------------------------------------------------------+
| ac/4+| OBSFRST|  Start-date of observations used in fit, YYYYMMDD         |
+------+--------+-----------------------------------------------------------+
| ac/4+| OBSLAST|  Stop-date of observations used in fit, YYYYMMDD          |
+------+--------+-----------------------------------------------------------+
| ac/5+|  PRELTV|  Planet relativity "bit-switch" byte: bits 0-7 are set to |
|      |        |  1 if relativity for corresponding planet was computed,   |
|      |        |  0 if not. For example, if Earth & Jupiter, FORTRAN(95)   |
|      |        |  statement IBITS(PRELTV,J,1) should return 1 when J=2 or  |
|      |        |  J=4, but zero for every other J through 7. No provision  |
|      |        |  for supporting Pluto relativity.                         |
+------+--------+-----------------------------------------------------------+
| ac/5+|  SPHMX3|  Earth grav. model max. degree; 0=point-mass, 2= J2 only, |
|      |        |  3= up to J3 zonal, 22= 2x2 field, 33=3x3 field, etc.     |
+------+--------+-----------------------------------------------------------+
| ac/5+|  SPHMX5|  Jupiter grav. max. deg.; 0=point-mass, 2= J2 only,       |
|      |        |  3= up to J3 zonal, 22= 2x2 field, 33=3x3 field, etc.     |
+------+--------+-----------------------------------------------------------+
| ac/5+|   JGSEP|  Galilean satellites used as sep. perturbers; 0=no 1=yes  |
+------+--------+-----------------------------------------------------------+
| ac/5+|  TWOBOD|  Two-body orbit model flag; 0=no 1=yes                    |
+------+--------+-----------------------------------------------------------+
| ac/5+|   NSATS|  Number of satellites; 99 if unknown.                     |
+------+--------+-----------------------------------------------------------+
| ac/4+|   UPARM|  Orbit condition code; 99 if not computed                 |
+------+--------+-----------------------------------------------------------+
| ac/4+|    LSRC|  Length of square-root cov. vector SRC (# elements used)  |
+------+--------+-----------------------------------------------------------+
| c/3+ |    IPYR|  Perihelion year (i.e., 1976, 2012, 2018, etc.)           |
+------+--------+-----------------------------------------------------------+
| ac/3+|    NDEL|  Number of radar delay measurements used in orbit soln.   |
+------+--------+-----------------------------------------------------------+
| ac/3+|    NDOP|  Number of radar Doppler measurements used in orbit soln. |
+------+--------+-----------------------------------------------------------+
| c/5+ |  NOBSMT|  Number of magnitude measurements used in total mag. soln.|
+------+--------+-----------------------------------------------------------+
| c/5+ |  NOBSMN|  Number of magnitude measurements used in nuc. mag. soln. |
+------+--------+-----------------------------------------------------------+
| c/3+ |  COMNUM|  IAU comet number (parsed from DESIG)                     |
+------+--------+-----------------------------------------------------------+
| ac/3+|  EQUNOX|  Equinox of orbital elements ('1950' or '2000')           |
+------+--------+-----------------------------------------------------------+
| ac/4+|   PENAM|  Planetary ephemeris ID/name                              |
+------+--------+-----------------------------------------------------------+
| ac/3+|   SBNAM|  Small-body perturber ephemeris ID/name                   |
+------+--------+-----------------------------------------------------------+
| a /3 |  SPTYPT|  Tholen spectral type                                     |
+------+--------+-----------------------------------------------------------+
| a /4+|  SPTYPS|  SMASS-II spectral type                                   |
+------+--------+-----------------------------------------------------------+
| ac/3+|    DARC|  Data arc span (year-year, OR integer # of days)          |
+------+--------+-----------------------------------------------------------+
| a /3+|  COMNT1|  Asteroid comment line #1                                 |
+------+--------+-----------------------------------------------------------+
| a /3+|  COMNT2|  Asteroid comment line #2                                 |
+------+--------+-----------------------------------------------------------+
| c/3+ |  COMNT3|  Comet comment line #1                                    |
+------+--------+-----------------------------------------------------------+
| c/3+ |  COMNT4|  Comet comment line #2                                    |
+------+--------+-----------------------------------------------------------+
| ac/3+|   DESIG|  Object designation                                       |
+------+--------+-----------------------------------------------------------+
| ac/4+|    ESTL|  Dynamic parameter estimation list. Last symbol set       |
|      |        |  to '+' if list is too long for field; check              |
|      |        |  object record comments field for full list.              |
+------+--------+-----------------------------------------------------------+
| ac/3+|    IREF|  Solution reference/ID/name                               |
+------+--------+-----------------------------------------------------------+
| ac/3+|    NAME|  Object name                                              |
+------+--------+-----------------------------------------------------------+