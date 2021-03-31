import erfa
from astropy import units as u
from astropy.coordinates.solar_system import PLAN94_BODY_NAME_TO_PLANET_INDEX

from ..constants import J2000
from ..frames import Planes
from .states import RVState


def get_mean_elements(body, epoch=J2000):
    """Get ecliptic mean elements of body.

    Parameters
    ----------
    body : ~poliastro.bodies.SolarSystemPlanet
        Body.
    epoch : astropy.time.Time
        Epoch.

    """
    # Internally, erfa.plan94 has a table of classical elements,
    # which are then converted to position and velocity...
    # So we are reverting the conversion in the last line
    # This way at least we avoid copy pasting the values
    try:
        name = body.name.lower()
        if name == "earth":
            name = "earth-moon-barycenter"

        body_index = PLAN94_BODY_NAME_TO_PLANET_INDEX[name]
        body_pv_helio = erfa.plan94(epoch.jd1, epoch.jd2, body_index)

        r = body_pv_helio["p"] * u.au
        v = body_pv_helio["v"] * u.au / u.day

    except KeyError as e:
        raise ValueError(f"No available mean elements for {body}") from e

    return RVState(body.parent, r, v, plane=Planes.EARTH_ECLIPTIC).to_classical()
