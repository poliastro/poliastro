import json

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from astropy.time import Time

from poliastro.contrib.czml_extract_default_params import DEFAULTS


class ExtractorCZML:
    """A class for extracting orbitary data to Cesium"""

    def __init__(self, orbit, start_epoch, end_epoch, N, fixed_rf=True):
        """
        Orbital constructor

        Parameters
        ----------
        orbit: poliastro.Orbit
            Orbit to be extracted
        start_epoch: ~astropy.time.core.Time
            Starting epoch
        end_epoch: ~astropy.time.core.Time
            Ending epoch
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
        self.start_epoch = ExtractorCZML.format_date(start_epoch.value)
        self.end_epoch = ExtractorCZML.format_date(end_epoch.value)
        self.init_czml()

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
            self.parse_dict_tuples([i] + t_key[1:], t_val)

        start_epoch = ExtractorCZML.format_date(
            min(Time(self.orbits[i][3]), Time(self.end_epoch)).iso
        )

        self.parse_dict_tuples(
            [i], [("id", str(i)), ("availability", start_epoch + "/" + self.end_epoch)]
        )
        self.parse_dict_tuples(
            [i, "path", "show"], [("interval", start_epoch + "/" + self.end_epoch)]
        )

        self.parse_dict_tuples(
            [i, "position"],
            [
                ("interpolationAlgorithm", "LAGRANGE"),
                ("interpolationDegree", 5),
                ("referenceFrame", "FIXED"),
                ("epoch", start_epoch),
                ("cartesian", list()),
            ],
        )
        self.init_orbit_packet_cords(i)

    def init_orbit_packet_cords(self, i):
        """

        Parameters
        ----------
        i: int
            Index of referenced orbit
        """
        h = (Time(self.end_epoch) - self.orbits[i][3]).to(u.second) / self.orbits[i][1]

        for k in range(self.orbits[i][1] + 2):
            cords = (
                self.orbits[i][0]
                .represent_as(CartesianRepresentation)
                .xyz.to(u.meter)
                .value
            )
            cords = np.insert(cords, 0, h.value * k, axis=0)

            self.czml[i]["position"]["cartesian"] += cords.tolist()
            self.orbits[i][0] = self.orbits[i][0].propagate(h)

    def init_czml(self):
        """
        Only called at the initialization of the extractor
        Builds packets.
        """

        self.parse_dict_tuples(
            [0], [("id", "document"), ("name", "simple"), ("version", "1.0")]
        )
        self.parse_dict_tuples(
            [0, "clock"],
            [
                ("interval", self.start_epoch + "/" + self.end_epoch),
                ("currentTime", self.start_epoch),
                ("multiplier", 60),
                ("range", "LOOP_STOP"),
                ("step", "SYSTEM_CLOCK_MULTIPLIER"),
            ],
        )

        self.init_orbit_packet(1)

    def change_id_params(self, i, id=None, name=None, description=None):
        """
        Change the id parameters.

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
            self.parse_dict_tuples([i], [("id", id)])
        if name is not None:
            self.parse_dict_tuples([i], [("name", name)])
        if description is not None:
            self.parse_dict_tuples([i], [("description", description)])

    def change_path_params(
        self, i, pixel_offset=None, color=None, width=None, show=None
    ):
        """
        Changes the path parameters.

        Parameters
        ----------
        i : int
            Referred body (count starts at 1)
        pixel_offset: list (int)
            The pixel offset (up and right)
        color: list (int)
            Rgba path color
        width: int
            Path width
        show: bool
            Indicates whether the path is visible
        """
        if pixel_offset is not None:
            self.parse_dict_tuples(
                [i, "label", "pixelOffset"], [("cartesian2", pixel_offset)]
            )
        if color is not None:
            self.parse_dict_tuples(
                [i, "path", "material", "solidColor", "color"], [("rgba", color)]
            )
        if width is not None:
            self.parse_dict_tuples([i, "path"], [("width", width)])
        if show is not None:
            self.parse_dict_tuples([i, "path", "show"], [("boolean", show)])

    def change_label_params(
        self, i, fill_color=None, outline_color=None, font=None, text=None, show=None
    ):
        """
        Change the label parameters.

        Parameters
        ----------
        i : int
            Referred body (count starts at 1)
        fill_color: list (int)
            Fill Color in rgba format
        outline_color: list (int)
            Outline Color in rgba format
        font: str
            Set label font style and size (CSS syntax)
        text: str
            Set label text
        show: bool
            Indicates whether the label is visible
        """
        if fill_color is not None:
            self.parse_dict_tuples([i, "label", "fillColor"], [("rgba", fill_color)])
        if outline_color is not None:
            self.parse_dict_tuples(
                [i, "label", "outlineColor"], [("rgba", outline_color)]
            )
        if font is not None:
            self.parse_dict_tuples([i, "label"], [("font", font)])
        if text is not None:
            self.parse_dict_tuples([i, "label"], [("text", text)])
        if show is not None:
            self.parse_dict_tuples([i, "label"], [("show", show)])

    def extract(self, ext_location=None):
        """
        Parameters
        ----------
        ext_location : str
            Path to extract your file to, if not path is given, return the dump in the console
        """

        return json.dumps(list(self.czml.values()))

    def add_orbit(self, orbit, N):
        """
        Adds an orbit

        Parameters
        ----------
        orbit: poliastro.Orbit
            Orbit to be added
        N: int
            Number of sample points
        """

        self.i += 1

        if orbit.epoch < Time(self.start_epoch):
            orbit = orbit.propagate(Time(self.start_epoch) - orbit.epoch)
        elif orbit.epoch > Time(self.end_epoch):
            raise ValueError(
                "The orbit's epoch cannot exceed the constructors ending epoch"
            )

        self.orbits[self.i] = [orbit, N, orbit.period, orbit.epoch]
        self.init_orbit_packet(self.i)

    def del_orbit(self, i):
        """
        Deletes an existing orbit

        Parameters
        ----------
        i: int
            Index of orbit to delete
        """

        del self.orbits[i]
        del self.czml[i]

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
        return date[:10] + "T" + date[11:-4] + "Z"
