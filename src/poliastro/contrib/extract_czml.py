import numpy as np
import json

from astropy.coordinates import CartesianRepresentation
from astropy import units as u

from poliastro.contrib.czml_extract_default_params import PIC_SATELLITE, DEFAULTS


class ExtractorCZML:
    """A class for extracting orbitary data to Cesium"""

    def __init__(self, orbit, N, fixed_rf=True):
        """
        Orbital constructor

        Parameters
        ----------
        orbit: poliastro.Orbit
            Orbit to be extracted
        N: int
            Number of sample points
        fixed_rf: bool
            Determines whether the reference frame is fixed
            When set to false, it assumes the frame moving
            with constant velocity.
        """
        self.czml = dict()
        self.orbits = {1: [orbit, N, orbit.period, orbit.epoch]}

        self.i = 1  # Current index id, used for insertion of new elements

        self.fixed_rf = fixed_rf
        self.start_epoch = ExtractorCZML.format_date(orbit.epoch.iso)
        self.end_epoch = ExtractorCZML.format_date((orbit.epoch + orbit.period).iso)
        self.init_czml()

        return None

    def parse_dict_tuples(self, path, tups):
        """
        Parameters
        ----------
        path : list (val)
            Dictionary path to insert to
        tups : list (val, val)
            Tuples to be assigned
        """

        # We only want to pass a reference czml, then by modifying our reference, we'll be modifying our base dictionary
        # which allows us to walk through a path of arbitrary length

        curr = self.czml
        for p in path:
            curr = curr.setdefault(p, {})
        for t in tups:
            curr[t[0]] = t[1]

    def init_orbit_packet(self, i):
        """
        Sets the default values for a given orbit

        Parameters
        ----------
        i : int
            Index of referenced orbit
        """
        for t_key, t_val in DEFAULTS:
            self.parse_dict_tuples(t_key, t_val)

        self.parse_dict_tuples(["b0", "clock"], [("interval", self.start_epoch + '/' + self.end_epoch),
                                                 ("currentTime", self.start_epoch), ("multiplier", 60)])
        self.parse_dict_tuples(["b" + str(i)], [("availability", self.start_epoch + '/' + self.end_epoch)])
        self.parse_dict_tuples(["b" + str(i), "path", "show"], [("interval", self.start_epoch + '/' + self.end_epoch)])

        self.parse_dict_tuples(["b" + str(i), "position"], [("interpolationAlgorithm", "LAGRANGE"),
                                                            ("interpolationDegree", 5), ("referenceFrame", "FIXED"),
                                                            ("epoch", self.start_epoch), ("cartesian", list())])
        self.init_orbit_packet_cords(i)

    def init_orbit_packet_cords(self, i):
        """

        Parameters
        ----------
        i: int
            Index of referenced orbit
        """
        h = self.orbits[i][2] / self.orbits[i][1]

        for k in range(self.orbits[i][1] + 1):
            cords = self.orbits[i][0].represent_as(CartesianRepresentation).xyz.to(u.meter).value
            cords = np.insert(cords, 0, h.value * k, axis=0)

            self.czml['b' + str(i)]['position']['cartesian'] += cords.tolist()
            self.orbits[i][0] = self.orbits[i][0].propagate(h)

    def init_czml(self):
        """
        Only called at the initialization of the extractor
        Builds packets.
        TODO: Parametrize the variables
        TODO: Add leadTime and trailTime
        """

        self.parse_dict_tuples(["b0"], [("id", "document"), ("name", "simple"), ("version", "1.0")])
        self.parse_dict_tuples(["b0", "clock"], [("interval", self.start_epoch + '/' + self.end_epoch),
                                                 ("currentTime", self.start_epoch), ("multiplier", 60)])

        self.init_orbit_packet(1)

    def change_id_params(self, i, id=None, name=None, description=None):
        """
        Parameters
        i : int
            Referred body (count starts at i)
        id: str
            Set orbit id
        name: str
            Set orbit name
        description: str
            Set orbit description
        """
        if id is not None:
            self.parse_dict_tuples(['b' + str(i)], [("id", id)])
        if name is not None:
            self.parse_dict_tuples(['b' + str(i)], [("name", name)])
        if description is not None:
            self.parse_dict_tuples(['b' + str(i)], [("description", description)])

    def change_label_params(self, i, fillColor=None, outlineColor=None, font=None, text=None, show=True):
        """
        Parameters
        n : int
            Referred body (count starts at 1)
        fillColor: list (int)
            Fill Color in rgba format
        outlineColor: list (int)
            Outline Color in rgba format
        font: str
            Set label font style and size (CSS syntax)
        text: str
            Set label text
        show: bool
            Indicates whether the label is visible
        """
        if fillColor is not None:
            self.parse_dict_tuples(['b' + str(i), "label", "fillColor"], [("rgba", fillColor)])
        if outlineColor is not None:
            self.parse_dict_tuples(['b' + str(i), "label", "outlineColor"], [("rgba", outlineColor)])
        if font is not None:
            self.parse_dict_tuples(['b' + str(i), "label"], [("font", font)])
        if text is not None:
            self.parse_dict_tuples(['b' + str(i), "label"], [("text", text)])
        if show is not None:
            self.parse_dict_tuples(['b' + str(i), "label"], [("show", show)])

    def extract(self, ext_location=None):
        """
        Parameters
        ----------
        ext_location : str
            Path to extract your file to, if not path is given, return the dump in the console
        """

        return json.dumps(list(self.czml.values()))

    @staticmethod
    def format_date(date):
        """
        Parameters
        ----------
        date : str
            date of the form "yyyy-mm-dd hh:mm:ss.ssss"

        Returns
        -------
        formatted_date : str
            date of the form "yyyy-mm-ddThh:mm:ssZ"
        """
        return date[:10] + 'T' + date[11:-4] + 'Z'
