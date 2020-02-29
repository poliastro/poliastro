"""
Catch that asteroid!
====================

"""


######################################################################
# First, we need to increase the timeout time to allow the download of
# data occur properly
# 

from astropy.utils.data import conf
conf.dataurl
print(conf.dataurl)

conf.remote_timeout 
print(conf.remote_timeout)

conf.remote_timeout = 10000


######################################################################
# Then, we do the rest of the imports and create our initial orbits.
# 

from astropy import units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris
import astropy.coordinates as coord
from astropy.coordinates import (
    GCRS,
    ICRS,
    CartesianDifferential,
    CartesianRepresentation,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
solar_system_ephemeris.set("jpl")

from poliastro.bodies import *
from poliastro.twobody import Orbit
from poliastro.plotting import StaticOrbitPlotter
from poliastro.plotting.misc import plot_solar_system

EPOCH = Time("2017-09-01 12:05:50", scale="tdb")

earth = Orbit.from_body_ephem(Earth, EPOCH)
print(earth)

earth.plot(label=Earth);
plt.show()

florence = Orbit.from_sbdb("Florence")
print(florence)


######################################################################
# Two problems: the epoch is not the one we desire, and the inclination is
# with respect to the ecliptic!
# 

print(florence.rv())

print(florence.epoch)

print(florence.epoch.iso)

print(florence.inc)


######################################################################
# We first propagate:
# 

florence = florence.propagate(EPOCH)
print(florence.epoch.tdb.iso)


######################################################################
# And now we have to convert to the same frame that the planetary
# ephemerides are using to make consistent comparisons, which is ICRS:
# 

def to_icrs(orbit):
    """Creates a new Orbit object with its coordinates transformed to ICRS.
    Notice that, strictly speaking, the center of ICRS is the Solar System Barycenter
    and not the Sun, and therefore these orbits cannot be propagated in the context
    of the two body problem. Therefore, this function exists merely for practical
    purposes.
    """

    coords = orbit.get_frame().realize_frame(orbit.represent_as(CartesianRepresentation, CartesianDifferential))
    coords.representation_type = CartesianRepresentation
    icrs_cart = coords.transform_to(ICRS).represent_as(CartesianRepresentation, CartesianDifferential)

    # Caution: the attractor is in fact the Solar System Barycenter
    ss = Orbit.from_vectors(
        Sun, r=icrs_cart.xyz, v=icrs_cart.differentials["s"].d_xyz, epoch=orbit.epoch
    )
    ss._frame = ICRS()
    return ss

florence_icrs = to_icrs(florence)
print(florence_icrs.rv())


######################################################################
# Let us compute the distance between Florence and the Earth:
# 

from poliastro.util import norm

print(norm(florence_icrs.r - earth.r) - Earth.R)


######################################################################
# .. raw:: html
# 
#    <div class="alert alert-success">
# 
# This value is consistent with what ESA says! :math:`7\,060\,160` km
# 
# .. raw:: html
# 
#    </div>
# 

abs(((norm(florence_icrs.r - earth.r) - Earth.R) - 7060160 * u.km) / (7060160 * u.km))

from IPython.display import HTML

HTML(
"""<blockquote class="twitter-tweet" data-lang="en"><p lang="es" dir="ltr">La <a href="https://twitter.com/esa_es">@esa_es</a> ha preparado un resumen del asteroide <a href="https://twitter.com/hashtag/Florence?src=hash">#Florence</a> 😍 <a href="https://t.co/Sk1lb7Kz0j">pic.twitter.com/Sk1lb7Kz0j</a></p>&mdash; AeroPython (@AeroPython) <a href="https://twitter.com/AeroPython/status/903197147914543105">August 31, 2017</a></blockquote>
<script src="//platform.twitter.com/widgets.js" charset="utf-8"></script>"""
)


######################################################################
# And now we can plot!
# 

frame = plot_solar_system(outer=False, epoch=EPOCH)
frame.plot(florence_icrs, label="Florence");
plt.show()


######################################################################
# The difference between doing it well and doing it wrong is clearly
# visible:
# 

frame = StaticOrbitPlotter()

frame.plot(earth, label="Earth")

frame.plot(florence, label="Florence (Ecliptic)")
frame.plot(florence_icrs, label="Florence (ICRS)");
plt.show()


######################################################################
# We can express Florence's orbit as viewed from Earth. In order to do
# that, we must set the Earth as the new attractor by making use of the
# ``change_attractor()`` method. However Florence is out of Earth's SOI,
# meaning that changing the attractor from Sun to Earth has no physical
# sense. We will make use of ``force=True`` argument so this method runs
# even if we know that we are out of new attractor's SOI.
# 

florence_hyper = florence.change_attractor(Earth, force=True)


######################################################################
# Previous warning was raised since Florence's orbit as seen from Earth is
# hyperbolic. Therefore if user wants to propagate this orbit along time,
# there will be some point at which the asteroid is out of Earth's
# influence (if not already).
# 


######################################################################
# We now retrieve the ephemerides of the Moon, which are given directly in
# GCRS:
# 

moon = Orbit.from_body_ephem(Moon, EPOCH)
print(moon)

moon.plot(label=Moon);
plt.show()


######################################################################
# And now for the final plot:
# 

import matplotlib.pyplot as plt

frame = StaticOrbitPlotter()

# This first plot sets the frame
frame.plot(florence_hyper, label="Florence")

# And then we add the Moon
frame.plot(moon, label=Moon)

plt.xlim(-1000000, 8000000)
plt.ylim(-5000000, 5000000)

plt.gcf().autofmt_xdate()
plt.show()


######################################################################
# .. raw:: html
# 
#    <div style="text-align: center; font-size: 3em;">
# 
# Per Python ad astra!
# 
# .. raw:: html
# 
#    </div>
# 
