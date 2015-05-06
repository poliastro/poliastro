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

SPK_TOP_URL = "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"
SPK_OLD_URL = "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/"

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


def select_kernel():
    """Selects appropriate kernel.

    Returns DE421 if available in data directory, else the first kernel found,
    else None.

    """
    kernel_files = glob.glob(os.path.join(SPK_LOCAL_DIR, "*.bsp"))
    if "de421.bsp" in kernel_files:
        kernel = SPK.open("de421.bsp")
    elif kernel_files:
        kernel = SPK.open(kernel_files[0])
    else:
        warnings.warn("""No SPICE kernels found under ~/.poliastro. \
Please download them manually or using

  poliastro download-spk [-d NAME]

to provide a default kernel, else pass a custom one as \
an argument to `planet_ephem`.""")
        kernel = None

    return kernel

default_kernel = select_kernel()


def download_kernel(name):
    """Downloads SPICE kernel by name.

    The function will try the top SPK path first, and then the old versions
    path in case the .bsp file is not found.

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
            print("Downloading %s.bsp from %s, please wait..." % (name, url_dir))
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
