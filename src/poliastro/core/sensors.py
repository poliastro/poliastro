import numpy as np
from numba import njit as jit


@jit
def min_and_max_ground_range(h, η_fov, η_center, R):
    """Calculates the minimum and maximum values of ground-range angles.

    Parameters
    ----------
    h : float
        Altitude over surface.
    η_fov : float
        Angle of the total area that a sensor can observe.
    η_center : float
        Center boresight angle.
    R : float
        Attractor equatorial radius.

    Returns
    -------
    Λ_min: float
        Minimum value of latitude and longitude.
    Λ_max: float
        Maximum value of latitude and longitude.

    Notes
    -----
    For further information, please take a look at "Fundamentals of Astrodynamics and Applications", 4th ed (2013)"
    by David A. Vallado, pages 853-860.

    """
    r_sat = R + h
    η_max = η_center + η_fov / 2
    η_min = η_center - η_fov / 2
    γ_max = np.arcsin(r_sat * np.sin(η_max) / R)

    if abs(γ_max) <= np.pi / 2:
        γ_max = np.pi - γ_max
    γ_min = np.arcsin(r_sat * np.sin(η_min) / R)
    if abs(γ_min) <= np.pi / 2:
        γ_min = np.pi - γ_min

    # Maximum and minimum slant ranges
    ρ_max = R * np.cos(γ_max) + r_sat * np.cos(η_max)
    ρ_min = R * np.cos(γ_min) + r_sat * np.cos(η_min)
    Λ_max = np.arcsin(ρ_max * np.sin(η_max) / R)
    Λ_min = np.arcsin(ρ_min * np.sin(η_min) / R)

    return Λ_min, Λ_max


@jit
def ground_range_diff_at_azimuth(h, η_fov, η_center, β, φ_nadir, λ_nadir, R):
    """Calculates the difference in ground-range angles from the η_center angle and the latitude and longitude of the target
    for a desired phase angle, β, used to specify where the sensor is looking.

    Parameters
    ----------
    h : float
        Altitude over surface.
    η_fov : float
        Angle of the total area that a sensor can observe.
    η_center : float
        Center boresight angle.
    β : float
        Phase angle, used to specify where the sensor is looking.
    φ_nadir : float
        Latitude angle of nadir point.
    λ_nadir : float
        Longitude angle of nadir point.
    R : float
        Earth equatorial radius.

    Returns
    -------
    delta_λ : float
        The difference in ground-range angles from the eta_center angle.
    φ_tgt: float
        Latitude angle of the target point.
    λ_tgt: float
        Longitude angle of the target point.

    Notes
    -----
    For further information, please take a look at "Fundamentals of Astrodynamics and Applications", 4th ed (2013)"
    by David A. Vallado, pages 853-860.

    """
    if not 0 <= β < np.pi:
        raise ValueError("beta must be between 0º and 180º")

    r_sat = R + h
    γ = np.arcsin(r_sat * np.sin(η_center) / R)
    if abs(γ) <= np.pi / 2:
        γ = np.pi - γ

    ρ = R * np.cos(γ) + r_sat * np.cos(η_center)
    Λ = np.arcsin(ρ * np.sin(η_center) / R)
    φ_tgt = np.arcsin(
        np.cos(β) * np.cos(φ_nadir) * np.sin(Λ) + np.sin(φ_nadir) * np.cos(Λ)
    )
    delta_Λ = np.arcsin(np.sin(β) * np.sin(Λ) / np.cos(φ_tgt))
    λ_tgt = λ_nadir + delta_Λ
    Λ_min, Λ_max = min_and_max_ground_range(h, η_fov, η_center, R)
    delta_λ = (Λ_max - Λ_min) / 2

    return delta_λ, φ_tgt, λ_tgt
