"""Function helpers.

"""
import numpy as np
from astropy import units as u
from astropy.time import Time
from numpy.linalg import norm as norm_np

from .core.util import alinspace as alinspace_fast


def norm(vec):
    """Norm of a Quantity vector that respects units.

    Parameters
    ----------
    vec : ~astropy.units.Quantity
        Vector with units.

    """
    return norm_np(vec.value) * vec.unit


def time_range(start, *, periods=50, spacing=None, end=None, format=None, scale=None):
    """Generates range of astronomical times.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    periods : int, optional
        Number of periods, default to 50.
    spacing : Time or Quantity, optional
        Spacing between periods, optional.
    end : Time or equivalent, optional
        End date.

    Returns
    -------
    Time
        Array of time values.

    """
    start = Time(start, format=format, scale=scale)

    if spacing is not None and end is None:
        result = start + spacing * np.arange(0, periods)

    elif end is not None and spacing is None:
        end = Time(end, format=format, scale=scale)
        result = start + (end - start) * np.linspace(0, 1, periods)

    else:
        raise ValueError("Either 'end' or 'spacing' must be specified")

    return result


@u.quantity_input(value=u.rad, values=u.rad)
def find_closest_value(value, values):
    """Calculates the closest value in the given values.

    Parameters
    ----------
    value : ~astropy.units.Quantity
        Reference value.
    values : ~astropy.units.Quantity
        Values to search from.

    """
    index = np.abs(np.asarray(values) * u.rad - value).argmin()
    return values[index]


@u.quantity_input(start=u.rad, stop=u.rad)
def alinspace(start, stop=None, *, num=50, endpoint=True):
    """Return increasing, evenly spaced angular values over a specified interval.

    """
    if stop is None:
        stop = start + 2 * np.pi * u.rad

    return (
        alinspace_fast(start.to_value(u.rad), stop.to_value(u.rad), num, endpoint)
        * u.rad
    )


def variable_step_bisection_algorithm(
    f,
    x0,
    y,
    step_size,
    args=(),
    kwargs={},
    xtol=None,
    xrtol=None,
    tol=None,
    rtol=None,
    maxiter=1000,
):

    r"""It calculates the closest value :math:`x` that approximate the function :math:`f` to the target value.
    This calculation process start with the initial value and tries to reduce the MSE error, :math:`(y - f(x))^2`,
    by adding a step size to :math:`x` in the direction that reduces the error.
    What makes this process interesting is that the step size is variable. How the step size changes is similar to
    VSLMS (Variable Sign Linear Mean Square) algorithm. The step size increases when the error has the same sign,
    avoiding multiple unnecessary iterations and when the error changes the sign it reduces the step size to obtain
    better approximation. The process will end when the error meets the tolerance conditions or when it exceeds the
    number of iterations.

    Parameters
    ----------
    f : function
        The function whose want to approximate to :math`y` value.
        It must be a function of a single variable of the form :math: `f(x,a,b,c...)`, where a,b,c...
        are extra arguments that can be passed in the args and kwargs parameters.
        It will use the form :math `f(x, *args, **kwargs)`
    x0 : float
        An initial estimate of :math:`x`. It's recommended to be near the actual solution to avoid extra computing
    y : float
        The value that want to be estimated
    step_size : float
        The initial value of the stepsize be added to :math:`x` in order to reduce the error :math:`|y - f(x)|`
    args : tuple, optional
        Extra arguments to be used in the function call.
    kwargs : dict, optional
        Extra arguments to be used in the function call.
    xtol : float, optional
        The allowable tolerance of the :math:`x`.
    xrtol : float, optional
        The allowable relative tolerance of the :math:`x`.
    tol : float, optional
        The allowable tolerance of the error, :math:`|y - f(x)|`.
    rtol : float, optional
        The allowable relative tolerance of the error using :math:`y` as reference.
    maxiter : int, optional
        Maximum number of iterations. Default value 1000.

    Returns
    -------
    x: float
        The closest obtained value of :math:`f` to :math:`y`. :math:`f(x) ~= y`

    Notes
    -----
    For further [information](https://iiav.org/ijav/content/volumes/21_2016_590031458046128/vol_1/835_fullpaper_1207561458214850.pdf)
    about the algorithm check the reference link.

    """
    # method used to determine if the end the algorithm
    def algorithm_ended(
        e, x, prev_x, y, step_size, itt, xtol, xrtol, tol, rtol, maxiter
    ):
        e, x, prev_x, y, step_size, xtol, tol = [
            k.value if isinstance(k, u.quantity.Quantity) else k
            for k in (e, x, prev_x, y, step_size, xtol, tol)
        ]
        # end because a solution has been met
        if e == 0:
            return True
        # end because step size is zero, there is nothing else to optimize
        if step_size == 0:
            return True
        # end because maximum precision has been reached
        if x == prev_x:
            return True
        # end by x tolerance
        if xtol is not None and step_size < xtol:
            return True
        # end by x relative tolerance
        if xrtol is not None and step_size < xrtol * abs(x):
            return True
        # end by error tolerance
        if tol is not None and e < tol:
            return True
        # end by relative error tolerance
        if rtol is not None and e < rtol * abs(y):
            return True
        # end by max iteration
        if maxiter is not None and itt > maxiter:
            return True
        # process must continue
        return False

    # calculate the estimation for x0 and find the stepsize direction where the error is reduced
    itt = 0
    prev_e = y - f(x0, *args, **kwargs)
    e = y - f(x0 + step_size, *args, **kwargs)
    step_size = (
        step_size
        if np.sign(prev_e) != np.sign(e) or abs(prev_e) > abs(e)
        else -step_size
    )
    # start iterative process
    prev_x = None
    x = x0
    e = prev_e
    e = e.value if isinstance(e, u.quantity.Quantity) else e
    step_increase = 1.0
    while not algorithm_ended(
        abs(e), x, prev_x, y, step_size, itt, xtol, xrtol, tol, rtol, maxiter
    ):
        prev_x = x.value if isinstance(x, u.quantity.Quantity) else x
        x += step_size
        prev_e = e
        e = y - f(x, *args, **kwargs)
        e = e.value if isinstance(e, u.quantity.Quantity) else e
        step_size = (
            step_size * step_increase
            if np.sign(e) == np.sign(prev_e)
            else -step_size / 2
        )
        # step increases only when there is not a previous sign change, to avoid a pivot state
        step_increase = 1.5 if np.sign(e) == np.sign(prev_e) else 1
        itt += 1
    return x


def calculate_ecc_from_a(a, nsol, inc, attractor):
    r"""Calculates the eccentricity for a ground-track orbit.

    Parameters
    ---------
    a: ~astropy.units.Quantity
        Semi-major axis.
    nsol: ~astropy.units.Quantity
        Mean rotation rate of the attractor around the sun.
    inc: ~astropy.units.Quantity
        Inclination.
    attractor: Body
        Main attractor.

    Returns
    -------
    ecc: ~astropy.units.Quantity
        Eccentricity.

    Notes
    -----

    .. math::
        F = \frac {-3\mu^{1/2}J_2R^{2}cos(inc)} {2a^{7/2}n_sol}
    .. math::
        ecc = (1 - (F((5sin^{2}(inc) - 4)cos(inc) + 2))^{1/2})^{1/2}

    where

    .. math::
        \partial{W}cos(i) + \partial{Ω} = n_{sol}

    Equation (11) from the paper; "Daily repeat-groundtrack Mars orbits", by Noreen, G. and Kerridge,
    S. and Diehl, R. and Neelon, J. and Ely, Todd and Turner, A.E. For further information please consult
    "Advances in the Astronautical Sciences", Vol 114, pages 1143-1155.

    """
    # Check if a is at least equal or greater than the minimum a using
    # dWa = dW*cos(inc) + dΩ = nsol
    # aux = -(3/4) * J2 * R**2 * µ**0.5
    # dWa = (aux * (5 * sin(inc) ** 2 - 4) * cos(inc) + 2 * aux * cos(inc)) / (a**3.5 * (1-ecc**2)**2) = nsol
    # (aux * (5 * sin(inc) ** 2 - 4) * cos(i) + 2 * aux * cos(inc)) / (nsol * (1-ecc**2)**2) = a**3.5
    # a >= (aux * (5 * sin(inc) ** 2 - 4) * cos(inc) + 2 * aux * cos(inc)) / nsol

    µ = attractor.k
    nsol = nsol.to_value(u.rad / u.s) / u.s
    aux = -(3 / 4) * attractor.J2 * (attractor.R ** 2) * (µ ** 0.5) * np.cos(inc)
    inc_factor = aux * (5 * (np.sin(inc) ** 2) - 4) + 2 * aux
    if (inc_factor / ((a ** (7 / 2)) * nsol)) < 0:
        ecc = (1 - 1e-10) * u.one
        a = (inc_factor / (nsol * (1 - ecc ** 2) ** 2)) ** (2 / 7)
    elif (inc_factor / ((a ** (7 / 2)) * nsol)) > 1:
        ecc = 0 * u.one
        a = (inc_factor / nsol) ** (2 / 7)
    else:
        ecc = (1 - (inc_factor / ((a ** (7 / 2)) * nsol)) ** 0.5) ** 0.5
    return a, ecc


def calculate_inc_from_a(a, nsol, ecc, attractor):

    r"""Calculates the inclination for a ground-track orbit.

    Parameters
    ---------
    a: ~astropy.units.Quantity
        Semi-major axis.
    nsol: ~astropy.units.Quantity
        Mean rotation rate of the attractor around the sun.
    ecc: ~astropy.units.Quantity
        Eccentricity.
    attractor: Body
        Main attractor.

    Returns
    -------
    inc: ~astropy.units.Quantity
        Inclination.

    Notes
    -----

    .. math::
        inc = cos^{-1}(\frac {2a^{7/2}n_{sol}(1 - ecc^{2})^{2}} {3J_{2}\sqrt\mu R^{2}})

    where

    .. math::
        \partial{Ω} = n_{sol}

    Equation (6) from the paper; "Daily repeat-groundtrack Mars orbits", by Noreen, G. and Kerridge,
    S. and Diehl, R. and Neelon, J. and Ely, Todd and Turner, A.E. For further information please consult
    "Advances in the Astronautical Sciences", Vol 114, pages 1143-1155.

    """
    # solve inc(a) by using approximation dΩ = nsol
    µ = attractor.k
    n = (µ / (a ** 3)) ** 0.5
    aux = (
        -0.75
        * n
        * attractor.J2.value
        * ((attractor.R.to(u.km) / (a * (1 - ecc ** 2))) ** 2)
    )
    dΩ = nsol
    cos_i = (dΩ / (2 * aux)).to_value(u.rad)
    cos_i = cos_i if abs(cos_i) <= 1 else cos_i / abs(cos_i)
    inc = np.arccos(cos_i) * u.rad
    return inc


def calculate_nt(a, nsol, attractor, ecc=None, inc=None):
    r"""Calculates the total mean motion of the orbit.

    Parameters
    ---------
    a: ~astropy.units.Quantity
        Semi-major axis.
    nsol: ~astropy.units.Quantity
        Mean rotation rate of the attractor around the sun.
    attractor: Body
        Main attractor.
    ecc: ~astropy.units.Quantity
        Eccentricity.
    inc: ~astropy.units.Quantity
        Inclination.


    Returns
    -------
    nt: ~astropy.units.Quantity
        Total mean motion

    Notes
    -----

    .. math::
        n_{t} = n + [\dot\Omega + (\dot\omega + \dot M_{o}) cos(i)] cos(i)

    Equation (5) from the paper; "Daily repeat-groundtrack Mars orbits", by Noreen, G. and Kerridge,
    S. and Diehl, R. and Neelon, J. and Ely, Todd and Turner, A.E. For further information please consult
    "Advances in the Astronautical Sciences", Vol 114, pages 1143-1155.

    """

    if ecc is None and inc is None:
        raise ValueError("Either ecc or inc or both values are needed")
    dΩ = None
    dW = None
    # specified inc but not ecc
    if ecc is None:
        _, ecc = calculate_ecc_from_a(a, nsol, inc, attractor)
        dΩ = nsol.to_value(u.rad / u.s) / u.s
        dW = 0
    # specified ecc but not inc
    µ = attractor.k
    n = (µ / (a ** 3)) ** 0.5
    aux = -0.75 * n * attractor.J2 * ((attractor.R / (a * (1 - ecc ** 2))) ** 2)
    if inc is None:
        # solve inc(a) by using approximation dΩ = nsol
        dΩ = nsol.to_value(u.rad / u.s) / u.s
        cos_i = dΩ / (2 * aux)
    else:
        cos_i = np.cos(inc.to_value(u.rad))
        dΩ = 2 * aux * cos_i
    sin2_i = 1 - cos_i ** 2
    dW = aux * (5 * sin2_i - 4) if dW is None else dW
    dM = aux * ((1 - ecc ** 2) ** 0.5) * (3 * sin2_i - 2)
    nt = n + (dΩ + (dW + dM) * cos_i) * cos_i
    return nt


def calculate_pq(a, nsol, attractor, ecc, inc):
    r"""Calculates the time it takes for an orbit to complete Q periods.

    Parameters
    ---------
    a: ~astropy.units.Quantity
        Semi-major axis.
    nsol: ~astropy.units.Quantity
        Mean rotation rate of the attractor around the sun.
    attractor: Body
        Main attractor.
    ecc: ~astropy.units.Quantity
        Eccentricity.
    inc: ~astropy.units.Quantity
        Inclination.

    Returns
    -------
    pq: ~astropy.units.Quantity
        Time it takes for an orbit to complete Q periods.

    Notes
    -----

    .. math::

        P_{q} = \frac {2 \pi}{n_{t}}

    Equation (10) from the paper; "Daily repeat-groundtrack Mars orbits", by Noreen, G. and Kerridge,
    S. and Diehl, R. and Neelon, J. and Ely, Todd and Turner, A.E. For further information please consult
    "Advances in the Astronautical Sciences", Vol 114, pages 1143-1155.

    """
    return (
        2 * np.pi / calculate_nt(a=a, nsol=nsol, attractor=attractor, ecc=ecc, inc=inc)
    )


def scipy_root_scalar(a, nsol, attractor, ecc, inc, pq_target):
    r""" Calculates the error between :math:`Pq` and :math:`Pq` approximation.

    Parameters:
    ----------
    a: ~astropy.units.Quantity
        Semi-major axis.
    nsol: ~astropy.units.Quantity
        Mean rotation rate of the attractor around the sun.
    attractor: Body
        Main attractor.
    ecc: ~astropy.units.Quantity
        Eccentricity.
    inc: ~astropy.units.Quantity
        Inclination.
    pq_target: ~astropy.units.Quantity
        Period to complete Q obits in a day.

    Returns
    ------
    result: ~astropy.units.Quantity
        Error.

    Notes
    -----

    .. math::

        P_{q} = \frac {2 \pi}{n_{t}}

    Equation (10) from the paper; "Daily repeat-groundtrack Mars orbits", by Noreen, G. and Kerridge,
    S. and Diehl, R. and Neelon, J. and Ely, Todd and Turner, A.E. For further information please consult
    "Advances in the Astronautical Sciences", Vol 114, pages 1143-1155.

    """
    a = a if isinstance(a, u.quantity.Quantity) else a * u.km
    pq = calculate_pq(a, nsol, attractor, ecc, inc)
    result = pq_target - pq
    return result.to_value(u.s)


def calculate_a_brackets(attractor, inc=None, ecc=None, Psol=None, Q=None, nsol=None):
    r""" Returns the bracket where there is a change of sign of the function in order to apply a root finding algorithm

    Parameters:
    -----------
    attractor: Body
        Main attractor.
    ecc: ~astropy.units.Quantity
        Eccentricity
    Psol: ~astropy.units.Quantity
        Solar Period of the attractor.
    Q: float
        Trace repetition parameter for a ground-track orbit,amount of orbits in a day, revolutions per day.
    nsol: ~astropy.units.Quantity
        Mean rotation rate of the attractor around the sun.

    Returns
    -------
    range_a ~astropy.units.Quantity
        An interval

    """
    nt = 2 * np.pi * Q / Psol
    min_a = attractor.R
    a = min_a
    nt_aprox = calculate_nt(a, nsol, attractor, ecc=ecc, inc=inc)
    sign = nt > nt_aprox
    max_a = None
    while a < 1000 * attractor.R:
        a = a + attractor.R
        nt_aprox = calculate_nt(a, nsol, attractor, ecc=ecc, inc=inc)
        new_sign = nt > nt_aprox
        if new_sign != sign:
            max_a = a
            break
        min_a = a
    if max_a is None:
        raise ValueError("Bracket couldn't be found")
    return min_a.to_value(u.km), max_a.to_value(u.km)
