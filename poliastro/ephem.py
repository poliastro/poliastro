# coding: utf-8
"""Planetary ephemerides.

Ideas
-----
* Don't download any data by default
* Provide a CLI utility to download data, and raise exception if none available

Example::

    $ poliastro ephem  # Downloads default
    $ poliastro ephem --name de421  # Downloads specific kernel

* The kernels are downloaded to ~/.poliastro, using os.path.expanduser
* Create a way to assign module-level kernel
* Python modules are singletons http://stackoverflow.com/q/6255050
* Use it for planet_ephem, or provide an extra argument
* Use entry points for provide a function

"""
from astropy import units as u
import glob
import warnings

from jplephem.spk import SPK


# TODO: Convert to enum, introduce requirement for Python 2
# https://caremad.io/2013/07/setup-vs-requirement/
# https://pip.pypa.io/en/latest/user_guide.html#requirements-files
MERCURY = 1
VENUS = 2
EARTH = 3
MARS = 4
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
PLUTO = 9
SUN = 10


def select_kernel():
    """Selects appropriate kernel.

    Returns DE421 if available in data directory, else the first kernel found,
    else None.

    """
    kernel_files = glob.glob("*.bsp")
    if "de421.bsp" in kernel_files:
        kernel = SPK.open("de421.bsp")
    elif kernel_files:
        kernel = SPK.open(kernel_files[0])
    else:
        warnings.warn("No SPICE kernels found under ~/.poliastro. Please "
                      "download them manually or using "
                      "`poliastro ephem --name KERNEL_NAME` to provide a "
                      "default kernel, else pass a custom one as an argument "
                      "to `planet_ephem`.")
        kernel = None

    return kernel

default_kernel = select_kernel()


def planet_ephem(body, epoch, kernel=default_kernel):
    """Position and velocity vectors of a given planet at a certain time.

    The vectors are computed with respect to the Solar System barycenter.

    Parameters
    ----------
    body : int
        Planetary body.
    epoch : astropy.Time
        Computation time.

    Returns
    -------
    r, v : Quantity
        Position and velocity vectors.

    """
    r, v = kernel[0, body].compute_and_differentiate(epoch.jd1, epoch.jd2)
    r *= u.km
    v *= u.km / u.s
    return r, v
