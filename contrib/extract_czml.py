import numpy as np
import json

from astropy.coordinates import CartesianRepresentation
from astropy import units as u

PIC_SATELLITE = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVDhPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47WoH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9iaqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYII="

DEFAULTS = [
    (["b0"], [("id", "document"), ("name", "simple"), ("version", "1.0")]),
    (["b0", "clock"], [("interval", "0000-00-00T00:00:00Z/0000-00-00T00:00:00Z"),
                ("currentTime", "0000-00-00T00:00:00Z/0000-00-00T00:00:00Z"),
                ("multiplier", 60), ("range", "LOOP_STOP"), ("step", "SYSTEM_CLOCK_MULTIPLIER")]),
    (["b1"], [("id", "ID"), ("name", "NAME"), ("availability", "0000-00-00T00:00:00Z/0000-00-00T00:00:00Z"),
                ("description", "DESCRIPTION")]),
    (["b1", "billboard", "eyeOffset"], [("cartesian", [0, 0, 0])]),
    (["b1", "billboard"], [("horizontalOrigin", "CENTER"), ("verticalOrigin", "CENTER"), ("image", PIC_SATELLITE),
                ("scale", 1), ("show", True)]),
    (["b1", "billboard", "pixelOffset"], [("cartesian2", [0, 0])]),
    (["b1", "label", "fillColor"], [("rgba", [255, 255, 0, 255])]),
    (["b1", "label"], [("font", "11pt Lucida Console"), ("horizontalOrigin", "LEFT"), ("outlineWidth", 2),
                ("show", True), ("style", "FILL_AND_OUTLINE"), ("text", "Text"), ("verticalOrigin", "CENTER")]),
    (["b1", "label", "pixelOffset"], [("cartesian2", [6, -4])]),
    (["b1", "label", "outlineColor"], [("rgba", [2, 100, 0, 255])]),
    (["b1", "path", "show"], [("interval", "0000-00-00T00:00:00Z/0000-00-00T00:00:00Z"), ("boolean", True)]),
    (["b1", "path"], [("width", 1), ("resolution", 120)]),
    (["b1", "path", "material", "solidColor", "color"], [("rgba", [255, 255, 0, 255])]),
]


class ExtractorCZML:
    """A class for extracting orbitary data to Cesium"""

    def __init__(self, orbit, N):
        """
        Orbital constructor

        Parameters
        ----------
        orbit: poliastro.Orbit
            Orbit to be extracted
        N: int
            Number of sample points
        """
        self.czml = dict()
        self.orbit = orbit
        self.N = N

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

    def init_czml(self):
        """
        Builds packets.
        TODO: Parametrize the variables
        TODO: Add leadTime and trailTime
        """

        # TODO: Read from file instead of hardcoded list
        for t_key, t_val in DEFAULTS:
            self.parse_dict_tuples(t_key, t_val)

        self.parse_dict_tuples(["b0", "clock"], [("interval", self.start_epoch + '/' + self.end_epoch),
                                ("currentTime", self.start_epoch), ("multiplier", 60)])
        self.parse_dict_tuples(["b1"], [("availability", self.start_epoch + '/' + self.end_epoch)])
        self.parse_dict_tuples(["b1", "path", "show"], [("interval", self.start_epoch + '/' + self.end_epoch)])

        self.parse_dict_tuples(["b1", "position"], [("interpolationAlgorithm", "LAGRANGE"), ("interpolationDegree", 5),
                                ("referenceFrame", "INERTIAL"), ("epoch", self.start_epoch), ("cartesian", list())])

        h = self.orbit.period/self.N

        for i in range(self.N+2):
            cords = self.orbit.represent_as(CartesianRepresentation).xyz.to(u.meter).value
            cords = np.insert(cords, 0, h.value*i, axis=0)

            self.czml['b1']['position']['cartesian'] += cords.tolist()
            self.orbit = self.orbit.propagate(h)

    def change_id_params(self, n, id=None, name=None, description=None):
        """
        Parameters
        n : int
            Referred body (count starts at 0)
        id: str
            Set orbit id
        name: str
            Set orbit name
        description: str
            Set orbit description
        """
        if id is not None:
            self.parse_dict_tuples(['b'+str(n)], [("id", id)])
        if name is not None:
            self.parse_dict_tuples(['b' + str(n)], [("name", name)])
        if description is not None:
            self.parse_dict_tuples(['b' + str(n)], [("description", description)])

    def change_label_params(self, n, fillColor=None, outlineColor=None, font=None, text=None, show=True):
        """
        Parameters
        n : int
            Referred body (count starts at 0)
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
            self.parse_dict_tuples(['b' + str(n), "label", "fillColor"], [("rgba", fillColor)])
        if outlineColor is not None:
            self.parse_dict_tuples(['b' + str(n), "label", "outlineColor"], [("rgba", outlineColor)])
        if font is not None:
            self.parse_dict_tuples(['b' + str(n), "label"], [("font", font)])
        if text is not None:
            self.parse_dict_tuples(['b' + str(n), "label"], [("text", text)])
        if show is not None:
            self.parse_dict_tuples(['b' + str(n), "label"], [("show", show)])


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

