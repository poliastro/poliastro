from astropy import _erfa as erfa, units as u
from astropy.coordinates.solar_system import PLAN94_BODY_NAME_TO_PLANET_INDEX

from ..constants import J2000
from ..frames import Planes
from .states import RVState


def get_mean_elements(body, epoch=J2000):
    """Get mean elements of body.

    """
    # Internally, erfa.plan94 has a table of classical elements,
    # which are then converted to position and velocity...
    # So we are reverting the conversion in the last line
    # This way at least we avoid copy pasting the values
    name = body.name.lower()
    if name in ("earth", "moon"):
        name = "earth-moon-barycenter"

    body_index = PLAN94_BODY_NAME_TO_PLANET_INDEX[name]
    body_pv_helio = erfa.plan94(epoch.jd1, epoch.jd2, body_index)

    return RVState(
        body.parent,
        body_pv_helio["p"] * u.au,
        body_pv_helio["v"] * u.au / u.day,
        plane=Planes.EARTH_ECLIPTIC,
    ).to_classical()
