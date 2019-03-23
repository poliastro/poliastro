import copy
import json

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from astropy.time import Time, TimeDelta

from poliastro.contrib.czml_extract_default_params import DEFAULTS
from poliastro.twobody.propagation import propagate


class CZMLExtractor:
    """A class for extracting orbitary data to Cesium"""

    def __init__(self, start_epoch, end_epoch, N):
        """
        Orbital constructor

        Parameters
        ----------
        start_epoch: ~astropy.time.core.Time
            Starting epoch
        end_epoch: ~astropy.time.core.Time
            Ending epoch
        N: int
            Default number of sample points.
            Unless otherwise specified, the number
            of sampled data points will be N when calling
            add_orbit()
        """
        self.czml = dict()
        self.orbits = []
        self.N = N
        self.i = 0

        self.start_epoch = CZMLExtractor.format_date(start_epoch)
        self.end_epoch = CZMLExtractor.format_date(end_epoch)

        self._init_czml_()

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

    def _init_orbit_packet_(self, i):
        """
        Sets the default values for a given orbit

        Parameters
        ----------
        i : int
            Index of referenced orbit
        """

        self.czml[i] = copy.deepcopy(DEFAULTS)

        start_epoch = CZMLExtractor.format_date(
            min(self.orbits[i][2], self.start_epoch)
        )

        self.parse_dict_tuples(
            [i],
            [
                ("id", str(i)),
                ("availability", start_epoch.value + "/" + self.end_epoch.value),
            ],
        )
        self.parse_dict_tuples(
            [i, "path", "show"],
            [("interval", start_epoch.value + "/" + self.end_epoch.value)],
        )

        self.parse_dict_tuples(
            [i, "position"],
            [
                ("interpolationAlgorithm", "LAGRANGE"),
                ("interpolationDegree", 5),
                ("referenceFrame", "FIXED"),
                ("epoch", start_epoch.value),
                ("cartesian", list()),
            ],
        )
        self._init_orbit_packet_cords_(i)

    def _init_orbit_packet_cords_(self, i):
        """

        Parameters
        ----------
        i: int
            Index of referenced orbit
        """
        h = (self.end_epoch - self.orbits[i][2]).to(u.second) / self.orbits[i][1]

        for k in range(self.orbits[i][1] + 2):
            position = propagate(self.orbits[i][0], TimeDelta(k * h))

            cords = position.represent_as(CartesianRepresentation).xyz.to(u.meter).value
            cords = np.insert(cords, 0, h.value * k, axis=0)

            self.czml[i]["position"]["cartesian"] += list(
                map(lambda x: x[0], cords.tolist())
            )

    def _init_czml_(self):
        """
        Only called at the initialization of the extractor
        Builds packets.
        """

        self.parse_dict_tuples(
            [-1], [("id", "document"), ("name", "simple"), ("version", "1.0")]
        )
        self.parse_dict_tuples(
            [-1, "clock"],
            [
                ("interval", self.start_epoch.value + "/" + self.end_epoch.value),
                ("currentTime", self.start_epoch.value),
                ("multiplier", 60),
                ("range", "LOOP_STOP"),
                ("step", "SYSTEM_CLOCK_MULTIPLIER"),
            ],
        )

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
        if ext_location:
            with open(ext_location, "w+") as fp:
                fp.write(json.dumps(list(self.czml.values())))

        return json.dumps(list(self.czml.values()))

    def add_orbit(self, orbit, N=None):
        """
        Adds an orbit

        Parameters
        ----------
        orbit: poliastro.Orbit
            Orbit to be added
        N: int
            Number of sample points
        """

        if N is None:
            N = self.N

        if orbit.epoch < Time(self.start_epoch):
            orbit = orbit.propagate(self.start_epoch - orbit.epoch)
        elif orbit.epoch > Time(self.end_epoch):
            raise ValueError(
                "The orbit's epoch cannot exceed the constructor's ending epoch"
            )

        self.orbits.append([orbit, N, orbit.epoch])

        self._init_orbit_packet_(self.i)

        self.i += 1

    @staticmethod
    def format_date(date):
        """
        Parameters
        ----------
        date : ~astropy.time.core.Time
            input date

        Returns
        -------
        formatted_date : ~astropy.time.core.Time
            ISO 8601 - compliant date
        """
        return Time(date, format="isot")
