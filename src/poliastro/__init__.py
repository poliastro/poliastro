"""
=========
poliastro
=========

Utilities and Python wrappers for Orbital Mechanics

"""

__version__ = "0.14.post0.dev0"

from .atmosphere import COESA62, COESA76
from .twobody import (Orbit, cowell, propagate, state_from_vector, D_to_nu, nu_to_D, nu_to_E, nu_to_F, E_to_nu,
                      F_to_nu, M_to_E, M_to_F, M_to_D, E_to_M, F_to_M, D_to_M, M_to_nu, nu_to_M, fp_angle,
                      raan_from_ltan, BaseState, ModifiedEquinoctialState, RVState, ClassicalState)
from .frames import (Planes, HeliocentricEclipticJ2000, GeocentricSolarEcliptic, GeocentricMeanEcliptic,
                     gcrs_to_geosolarecliptic, geosolarecliptic_to_gcrs, ICRS, HCRS, MercuryICRS, VenusICRS, GCRS,
                     MarsICRS, JupiterICRS, SaturnICRS, UranusICRS, NeptuneICRS, PlutoICRS, MoonICRS, SunFixed,
                     MercuryFixed, VenusFixed, ITRS, MarsFixed, JupiterFixed, SaturnFixed, UranusFixed, NeptuneFixed,
                     PlutoFixed, MoonFixed, get_frame)
from .czml import ellipsoidal_to_cartesian, intersection_ellipsoid_line, project_point_on_ellipsoid
from .iod import lambert
from .neos import (asteroid_db, comet_db, orbit_from_name, orbit_from_record, record_from_name,
                   string_record_from_name, read_headers, read_record, download_dastcom5, entire_db)
from .plotting import (OrbitPlotter2D, OrbitPlotter3D, StaticOrbitPlotter, plot_solar_system, porkchop, generate_circle,
                       generate_sphere, generate_label)
from .twobody.thrust import change_a_inc, change_argp, change_ecc_quasioptimal, change_inc_ecc
from astropy import units as u
from astropy.time import Time
from .bodies import (Sun, Moon, SolarSystemBody, Body, _q, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus,
                     Neptune, Pluto)
from .examples import iss, molniya, _a, _r_a, _r_p, churi, soyuz_gto
from .core import (J2_perturbation, J3_perturbation, atmospheric_drag, shadow_function, third_body,radiation_pressure,
                   rv_pqw, coe_rotation_matrix, coe2mee, coe2rv, rv2coe, mee2coe, hyp2f1b, M_parabolic,
                   M_parabolic_prime, newton, D_to_nu, nu_to_D, nu_to_E, nu_to_F, E_to_nu, F_to_nu, M_to_E, M_to_F,
                   M_to_D, E_to_M, F_to_M, D_to_M, M_to_nu, nu_to_M, fp_angle, vallado, izzo, dct2, idct2, rsmooth,
                   bisquare, c2, c3, circular_velocity, rotation_matrix, mean_motion, kepler, mikkola, markley,
                   pimienta, gooding, danby, func_twobody)
from .threebody import compute_flyby, lagrange_points_vec, laplace_radius, lagrange_points, hill_radius
from .maneuver import Maneuver
from .coordinates import inertial_body_centered_to_pqw
from .ephem import build_ephem_interpolant
from .spheroid_location import SpheroidLocation
from .constants import *
from .util import (circular_velocity, norm, time_range, hyp_nu_limit, get_eccentricity_critical_argp,
                   get_eccentricity_critical_inc, get_inclination_critical_argp, find_closest_value)

__all__ = [
    # .twobody
    "Orbit", "cowell", "propagate", "state_from_vector", "D_to_M", "D_to_nu", "E_to_M", "E_to_nu", "fp_angle",
    "F_to_M", "F_to_nu", "M_to_nu", "M_to_D", "M_to_F", "M_to_E", "raan_from_ltan", "nu_to_M", "nu_to_F", "nu_to_E",
    "nu_to_D", "BaseState", "ModifiedEquinoctialState", "RVState", "ClassicalState",

    # .atmosphere
    "COESA62", "COESA76",

    # .core
    "J2_perturbation", "J3_perturbation", "atmospheric_drag", "radiation_pressure", "shadow_function", "third_body", "rv_pqw", "coe2rv", "coe2mee", "coe_rotation_matrix", "mee2coe", "rv2coe",
    "M_parabolic", "M_parabolic_prime", "newton", "D_to_nu", "nu_to_D", "nu_to_E", "nu_to_F", "E_to_nu", "F_to_nu",
    "M_to_E", "M_to_F", "M_to_D", "E_to_M", "F_to_M", "D_to_M", "M_to_nu", "nu_to_M", "fp_angle", "hyp2f1b",
    "vallado", "izzo", "dct2", "idct2", "rsmooth", "bisquare", "c2", "c3", "circular_velocity", "rotation_matrix",
    "mean_motion", "kepler", "mikkola", "markley", "pimienta", "gooding", "danby", "func_twobody",



    # .czml
    "ellipsoidal_to_cartesian", "intersection_ellipsoid_line", "project_point_on_ellipsoid",

    # .frames
    "Planes", "HeliocentricEclipticJ2000", "GeocentricMeanEcliptic", "GeocentricSolarEcliptic",
    "gcrs_to_geosolarecliptic", "geosolarecliptic_to_gcrs", "ICRS", "HCRS", "MercuryICRS", "VenusICRS", "GCRS",
    "MarsICRS", "JupiterICRS", "SaturnICRS", "UranusICRS", "NeptuneICRS", "PlutoICRS", "MoonICRS",
    "SunFixed", "MercuryFixed", "VenusFixed", "ITRS", "MarsFixed", "JupiterFixed", "SaturnFixed", "UranusFixed",
    "NeptuneFixed", "PlutoFixed", "MoonFixed", "get_frame",

    # .iod
    "lambert",

    # .neos
    "asteroid_db", "comet_db", "orbit_from_record", "orbit_from_name", "read_record", "read_headers",
    "record_from_name", "entire_db", "string_record_from_name", "download_dastcom5",

    # .plotting
    "OrbitPlotter2D", "OrbitPlotter3D", "StaticOrbitPlotter", "plot_solar_system", "porkchop", "generate_label",
    "generate_sphere", "generate_circle",

    # .threebody
    "compute_flyby", "lagrange_points", "laplace_radius", "lagrange_points_vec", "hill_radius",

    # .twobody.thrust
    "change_a_inc", "change_argp", "change_ecc_quasioptimal", "change_inc_ecc",

    # from astropy
    "u", "Time",

    # .bodies
    "Sun", "Moon", "SolarSystemBody", "Mercury", "Venus", "Earth", "Mars", "Jupiter",
    "Saturn", "Uranus", "Neptune", "Pluto", "Body", "_q",

    # .examples
    "iss", "churi", "molniya", "_a", "_r_p", "_r_a", "soyuz_gto",

    # .maneuver
    "Maneuver",

    # .coordinates
    'inertial_body_centered_to_pqw',

    # .ephem
    "build_ephem_interpolant",

    # .spheroid_location
    "spheroid_location",

    # .utils
    "circular_velocity", "norm", "time_range", "hyp_nu_limit", "get_eccentricity_critical_argp",
    "get_eccentricity_critical_inc", "get_inclination_critical_argp", "find_closest_value",

    # .constants
    "J2000",
    "J2000_TDB",
    "J2000_TT",
    "GM_sun",
    "GM_earth",
    "GM_mercury",
    "GM_venus",
    "GM_mars",
    "GM_jupiter",
    "GM_saturn",
    "GM_uranus",
    "GM_neptune",
    "GM_pluto",
    "GM_moon",
    "M_earth",
    "M_jupiter",
    "M_sun",
    "R_mean_earth",
    "R_mean_mercury",
    "R_mean_venus",
    "R_mean_mars",
    "R_mean_jupiter",
    "R_mean_saturn",
    "R_mean_uranus",
    "R_mean_neptune",
    "R_mean_pluto",
    "R_mean_moon",
    "R_earth",
    "R_mercury",
    "R_venus",
    "R_mars",
    "R_jupiter",
    "R_saturn",
    "R_sun",
    "R_uranus",
    "R_neptune",
    "R_pluto",
    "R_moon",
    "R_polar_earth",
    "R_polar_mercury",
    "R_polar_venus",
    "R_polar_mars",
    "R_polar_jupiter",
    "R_polar_saturn",
    "R_polar_uranus",
    "R_polar_neptune",
    "R_polar_pluto",
    "R_polar_moon",
    "J2_sun",
    "J2_earth",
    "J3_earth",
    "J2_mars",
    "J3_mars",
    "J2_venus",
    "J3_venus",
    "H0_earth",
    "rho0_earth",
    "Wdivc_sun",
    "rotational_period_earth",
    "rotational_period_sun",
    "rotational_period_mercury",
    "rotational_period_venus",
    "rotational_period_moon",
    "rotational_period_mars",
    "rotational_period_jupiter",
    "rotational_period_saturn",
    "rotational_period_uranus",
    "rotational_period_neptune",
    "rotational_period_pluto",
]
