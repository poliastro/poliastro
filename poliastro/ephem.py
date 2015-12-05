# coding: utf-8
"""Planetary ephemerides.

"""
import os
import glob
import warnings
warnings.formatwarning = lambda msg, *_: str(msg) + '\n'

from astropy import units as u

from jplephem.spk import SPK

try:
    from urllib.request import urlretrieve
    from urllib.error import HTTPError
except ImportError:  # Python 2
    from urllib import urlretrieve
    from urllib2 import HTTPError

NAIF_BASE_URL = "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/"
SPK_TOP_URL = NAIF_BASE_URL + "planets/"
SPK_OLD_URL = NAIF_BASE_URL + "planets/a_old_versions/"

SPK_LOCAL_DIR = os.path.expanduser("~/.poliastro")


# Planets
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


def select_kernel(kernel_filenames):
    """Selects appropriate kernel filename from a list.

    Returns DE421 if present, else the first name found, else None.

    .. versionadded:: 0.3.0

    """
    if "de421.bsp" in kernel_filenames:
        kernel = "de421.bsp"
    elif kernel_filenames:
        kernel = kernel_filenames[0]
    else:
        kernel = None

    return kernel

kernel_fname = select_kernel(
    glob.glob(os.path.join(SPK_LOCAL_DIR, "*.bsp")))

if kernel_fname:
    default_kernel = SPK.open(kernel_fname)
else:
    warnings.warn("""No SPICE kernels found under ~/.poliastro.
Please download them manually or using

  poliastro download-spk [-d NAME]

to provide a default kernel, else pass a custom one as
an argument to `planet_ephem`.""")
    default_kernel = None


def download_kernel(name):
    """Downloads SPICE kernel by name.

    The function will try the top SPK path first, and then the old versions
    path in case the .bsp file is not found.

    .. versionadded:: 0.3.0

    """
    destination_file = os.path.join(SPK_LOCAL_DIR, name + ".bsp")
    if os.path.isfile(destination_file):
        # Don't download file again
        print("File %s.bsp already exists under %s" % (name, SPK_LOCAL_DIR))
        return

    # Create .poliastro directory if not existing
    # Simple solution, no race condition checking, see
    # http://stackoverflow.com/a/273227
    if not os.path.isdir(SPK_LOCAL_DIR):
        os.makedirs(SPK_LOCAL_DIR)

    # Actually download file
    for url_dir in SPK_TOP_URL, SPK_OLD_URL:
        try:
            bsp_url = url_dir + name + ".bsp"
            print("Downloading %s.bsp from %s, please wait..." %
                  (name, url_dir))
            urlretrieve(bsp_url, destination_file)
        except HTTPError as e:
            print(e.msg)
            continue
        else:
            break
    else:
        raise RuntimeError("File %s.bsp not found under top nor old versions "
                           "directories, please consider downloading it "
                           "manually" % name)


def planet_ephem(body, epoch, kernel=default_kernel):
    """Position and velocity vectors of a given planet at a certain time.

    The vectors are computed with respect to the Solar System barycenter.

    .. versionadded:: 0.3.0

    Parameters
    ----------
    body : int
        Planetary body.
    epoch : astropy.time.Time
        Computation time. Can be scalar or vector.
    kernel : jplephem.spk.SPK, optional
        jplephem SPK kernel to make the computation, if not given a default
        one will be used.

    Returns
    -------
    r, v : Quantity
        Position and velocity vectors.

    """
    r, v = kernel[0, body].compute_and_differentiate(epoch.jd1, epoch.jd2)
    r = r * u.km
    v = v * u.km / u.day
    return r, v
