from datetime import timezone
from typing import Any, List

import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
from astropy.time import Time, TimeDelta
from czml3.core import Document, Packet, Preamble
from czml3.enums import InterpolationAlgorithms, ReferenceFrames
from czml3.properties import (
    Billboard,
    Clock,
    Color,
    Label,
    Material,
    Path,
    Position,
    SolidColorMaterial,
)
from czml3.types import IntervalValue, TimeInterval
from erfa import gd2gce

from poliastro.bodies import Earth
from poliastro.core.czml_utils import (
    project_point_on_ellipsoid as project_point_on_ellipsoid_fast,
)

PIC_SATELLITE = (
    "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAX"
    "NSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVD"
    "hPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47W"
    "oH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9i"
    "aqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp"
    "1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYI"
    "I="
)
PIC_GROUNDSTATION = (
    "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAX"
    "NSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACvSURBVD"
    "hPrZDRDcMgDAU9GqN0lIzijw6SUbJJygUeNQgSqepJTyHG91LVVpwDdfxM3T9TSl1EXZvDwi"
    "i471fivK73cBFFQNTT/d2KoGpfGOpSIkhUpgUMxq9DFEsWv4IXhlyCnhBFnZcFEEuYqbiUlN"
    "wWgMTdrZ3JbQFoEVG53rd8ztG9aPJMnBUQf/VFraBJeWnLS0RfjbKyLJA8FkT5seDYS1Qwyv"
    "8t0B/5C2ZmH2/eTGNNBgMmAAAAAElFTkSuQmCC"
)


class CZMLExtractor:
    """A class for extracting orbitary data to Cesium"""

    def __init__(
        self,
        start_epoch,
        end_epoch,
        N,
        attractor=None,
        pr_map=None,
        scene3D=True,
    ):
        """Orbital constructor

        Parameters
        ----------
        start_epoch : ~astropy.time.core.Time
            Starting epoch
        end_epoch : ~astropy.time.core.Time
            Ending epoch
        N : int
            Default number of sample points. Unless otherwise specified, the
            number of sampled data points will be N when calling add_orbit()
        attractor : poliastro.bodies.Body
            Attractor of the orbits
        pr_map : str
            A URL to the projection of the defined ellipsoid (UV map).
        scene3D : bool
            Determines the scene mode. If set to true, the scene is set to 3D
            mode, otherwise it's the orthographic projection.

        """
        self.packets = []  # type: List[Packet]
        self.trajectories = []  # type: List[Any]
        self.attractor = attractor
        self.orbits = []  # type: List[Any]
        self.N = N
        self.i = 0
        self.gs_n = 0

        if not self.attractor:
            self.attractor = Earth
        elif (
            self.attractor.R is None
            or self.attractor.R == 0
            or self.attractor.R_polar is None
            or self.attractor.R_polar == 0
        ):
            raise ValueError(
                "Invalid ellipsoid of attractor.\n"
                + "Make sure your body has valid 'R' and 'R_polar' parameters"
            )

        ellipsoid = (
            self.attractor.R.to_value(u.m),
            self.attractor.R.to_value(u.m),
            self.attractor.R_polar.to_value(u.m),
        )

        self.cust_prop = [ellipsoid, pr_map, scene3D]

        self.start_epoch = Time(start_epoch, format="isot")
        self.end_epoch = Time(end_epoch, format="isot")

        self._init_czml_()

        self._change_custom_params(*self.cust_prop)

    def _init_orbit_packet_cords_(self, i, rtol):
        """
        Parameters
        ----------
        i : int
            Index of referenced orbit
        rtol : float
            Maximum relative error permitted

        Returns
        -------
        coordinate list
        """
        cart_cords = []  # type: List[float]

        h = (self.end_epoch - self.orbits[i][2]).to(u.s) / self.orbits[i][1]

        # Get rounding factor given the relative tolerance
        rf = 0
        while rtol < 1:
            rtol *= 10
            rf += 1

        for k in range(self.orbits[i][1] + 2):
            # FIXME: Unused rtol with default propagation method,
            # should we just get rid of it?
            cords = (
                self.orbits[i][0]
                .propagate(TimeDelta(k * h))
                .r.to_value(u.m)[..., None]
            )

            cords = np.insert(cords, 0, h.value * k, axis=0)

            # Flatten list
            cart_cords += list(map(lambda x: round(x[0], rf), cords.tolist()))

        return cart_cords

    def _init_groundtrack_packet_cords_(self, i, rtol):
        """
        Parameters
        ----------
        i : int
            Index of referenced orbit
        rtol : float
            Maximum relative error permitted

        Returns
        -------
        coordinate list
        """
        cart_cords = []  # type: List[float]

        h = (self.end_epoch - self.orbits[i][2]).to(u.s) / self.orbits[i][1]

        # Get rounding factor given the relative tolerance
        rf = 0
        while rtol < 1:
            rtol *= 10
            rf += 1

        ellipsoid = self.cust_prop[0]

        for k in range(self.orbits[i][1] + 2):
            # FIXME: Unused rtol with default propagation method,
            # should we just get rid of it?
            cords = (
                self.orbits[i][0]
                .propagate(TimeDelta(k * h))
                .r.to_value(u.m)[..., None]
            )

            cords = np.insert(cords, 0, h.value * k, axis=0)

            # Flatten list
            cords = list(map(lambda x: round(x[0], rf), cords.tolist()))
            t, p = cords[0], cords[1:]
            pr_p = project_point_on_ellipsoid_fast(
                p[0], p[1], p[2], ellipsoid[0], ellipsoid[1], ellipsoid[2]
            )
            # Add a small number to ensure that our point lies above the surface of the
            # ellipsoid. We do this because small losses in precision may cause the point
            # to lie slightly below the surface. An iterative method could be used instead
            # but the error margin is too small to be worth it.

            _cords = t, pr_p[0] + 0.1, pr_p[1] + 0.1, pr_p[2] + 0.1
            cart_cords += _cords

        return cart_cords

    def _init_czml_(self):
        """
        Only called at the initialization of the extractor Builds packets.
        """
        pckt = Preamble(
            name="document_packet",
            clock=IntervalValue(
                start=self.start_epoch,
                end=self.end_epoch,
                value=Clock(
                    currentTime=self.start_epoch.datetime.replace(
                        tzinfo=timezone.utc
                    ),
                    multiplier=60,
                ),
            ),
        )
        self.packets.append(pckt)

    def _change_custom_params(self, ellipsoid, pr_map, scene3D):
        """
        Change the custom properties package.

        Parameters
        ----------
        ellipsoid : list[int]
            Defines the attractor ellipsoid. The list must have three numbers
            representing the radii in the x, y and z axis
        pr_map : str
            A URL to the projection of the defined ellipsoid (UV map).
        scene3D : bool
            If set to true, the scene is set to 3D mode,
            otherwise it's the orthographic projection.

        """

        if pr_map is None:
            pr_map = (
                "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg",
            )
            self.cust_prop[1] = pr_map

        custom_props = {
            "custom_attractor": True,
            "ellipsoid": [{"array": ellipsoid}],
            "map_url": pr_map,
            "scene3D": scene3D,
        }

        pckt = Packet(id="custom_properties", properties=custom_props)

        self.packets.append(pckt)

    def add_ground_station(
        self,
        pos,
        id_description=None,
        label_fill_color=None,
        label_font=None,
        label_outline_color=None,
        label_text=None,
        label_show=True,
    ):
        """
        Adds a ground station

        Parameters
        ----------
        pos : list[~astropy.units.Quantity]
            Coordinates of ground station,
            list of geodetic latitude and longitude [lon, lat] (0 elevation)
        id_description : str
            Set ground station description
        label_fill_color : list[int]
            Fill Color in rgba format
        label_outline_color : list[int]
            Outline Color in rgba format
        label_font : str
            Set label font style and size (CSS syntax)
        label_text : str
            Set label text
        label_show : bool
            Indicates whether the label is visible
        """
        if (
            len(pos) == 2
            and isinstance(pos[0], u.quantity.Quantity)
            and isinstance(pos[0], u.quantity.Quantity)
        ):
            u0, v0 = pos

            if self.cust_prop[0]:
                a, b = (
                    self.cust_prop[0][0],
                    self.cust_prop[0][2],
                )  # Get semi-major and semi-minor axises
            else:
                a, b = Earth.R.to_value(u.m), Earth.R_polar.to_value(u.m)

            f = 1 - (b / a)  # Flattenning

            pos = list(gd2gce(a, f, u0.to_value(u.rad), v0.to_value(u.rad), 0))

        else:
            raise TypeError(
                "Invalid coordinates. Coordinates must be of the form [u, v] where u, v are astropy units"
            )

        pckt = Packet(
            id="GS" + str(self.gs_n),
            description=id_description,
            availability=TimeInterval(
                start=self.start_epoch, end=self.end_epoch
            ),
            position=Position(cartesian=pos),
            label=Label(
                show=label_show,
                text=label_text,
                font=label_font
                if label_font is not None
                else "11pt Lucida Console",
                fillColor=Color(rgba=label_fill_color)
                if label_fill_color is not None
                else None,
                outlineColor=Color(rgba=label_outline_color)
                if label_outline_color is not None
                else None,
            ),
            billboard=Billboard(image=PIC_GROUNDSTATION, show=True),
        )

        self.packets.append(pckt)
        self.gs_n += 1

    def add_orbit(
        self,
        orbit,
        rtol=1e-10,
        N=None,
        groundtrack_show=False,
        groundtrack_lead_time=None,
        groundtrack_trail_time=None,
        groundtrack_width=None,
        groundtrack_color=None,
        id_name=None,
        id_description=None,
        path_width=None,
        path_show=None,
        path_color=None,
        label_fill_color=None,
        label_outline_color=None,
        label_font=None,
        label_text=None,
        label_show=None,
    ):
        """
        Adds an orbit

        Parameters
        ----------
        orbit : poliastro.twobody.orbit.Orbit
            Orbit to be added
        rtol : float
            Maximum relative error permitted
        N : int
            Number of sample points
        groundtrack_show : bool
            If set to true, the groundtrack is
            displayed.
        groundtrack_lead_time : float
            The time the animation is ahead of the real-time groundtrack
        groundtrack_trail_time : float
            The time the animation is behind the real-time groundtrack
        groundtrack_width : int
            Groundtrack width
        groundtrack_color : list[int]
            Rgba groundtrack color. By default, it is set to the path color
        id_name : str
            Set orbit name
        id_description : str
            Set orbit description
        path_width : int
            Path width
        path_show : bool
            Indicates whether the path is visible
        path_color : list[int]
            Rgba path color
        label_fill_color : list[int]
            Fill Color in rgba format
        label_outline_color : list[int]
            Outline Color in rgba format
        label_font : str
            Set label font style and size (CSS syntax)
        label_text : str
            Set label text
        label_show : bool
            Indicates whether the label is visible
        """

        if N is None:
            N = self.N

        if orbit.epoch < self.start_epoch:
            orbit = orbit.propagate(self.start_epoch - orbit.epoch)
        elif orbit.epoch > self.end_epoch:
            raise ValueError(
                "The orbit's epoch cannot exceed the constructor's ending epoch"
            )

        if rtol <= 0 or rtol >= 1:
            raise ValueError(
                "The relative tolerance must be a value in the range (0, 1)"
            )

        self.orbits.append([orbit, N, orbit.epoch])
        cartesian_cords = self._init_orbit_packet_cords_(self.i, rtol=rtol)

        start_epoch = Time(
            min(self.orbits[self.i][2], self.start_epoch), format="isot"
        )

        pckt = Packet(
            id=self.i,
            name=id_name,
            description=id_description,
            availability=TimeInterval(
                start=self.start_epoch, end=self.end_epoch
            ),
            position=Position(
                interpolationDegree=5,
                interpolationAlgorithm=InterpolationAlgorithms.LAGRANGE,
                referenceFrame=ReferenceFrames.INERTIAL,
                cartesian=cartesian_cords,
                # Use explicit UTC timezone, rather than the default, which is a local timezone.
                epoch=start_epoch.datetime.replace(tzinfo=timezone.utc),
            ),
            path=Path(
                show=path_show,
                width=path_width,
                material=Material(
                    solidColor=SolidColorMaterial(
                        color=Color.from_list(path_color)
                    )
                )
                if path_color is not None
                else Material(
                    solidColor=SolidColorMaterial(
                        color=Color.from_list([255, 255, 0])
                    )
                ),
                resolution=120,
            ),
            label=Label(
                text=label_text,
                font=label_font
                if label_font is not None
                else "11pt Lucida Console",
                show=label_show,
                fillColor=Color(rgba=label_fill_color)
                if label_fill_color is not None
                else Color(rgba=[255, 255, 0, 255]),
                outlineColor=Color(rgba=label_outline_color)
                if label_outline_color is not None
                else Color(rgba=[255, 255, 0, 255]),
            ),
            billboard=Billboard(image=PIC_SATELLITE, show=True),
        )

        self.packets.append(pckt)

        if groundtrack_show:

            groundtrack_color = path_color

            groundtrack_cords = self._init_groundtrack_packet_cords_(
                self.i, rtol=rtol
            )
            pckt = Packet(
                id="groundtrack" + str(self.i),
                availability=TimeInterval(
                    start=self.start_epoch, end=self.end_epoch
                ),
                position=Position(
                    interpolationDegree=5,
                    interpolationAlgorithm=InterpolationAlgorithms.LAGRANGE,
                    referenceFrame=ReferenceFrames.INERTIAL,
                    cartesian=groundtrack_cords,
                    # Use explicit UTC timezone, rather than the default, which is a local timezone.
                    epoch=start_epoch.datetime.replace(tzinfo=timezone.utc),
                ),
                path=Path(
                    show=True,
                    material=Material(
                        solidColor=SolidColorMaterial(
                            color=Color(rgba=groundtrack_color)
                        )
                    )
                    if groundtrack_color is not None
                    else Material(
                        solidColor=SolidColorMaterial(
                            color=Color(rgba=[255, 255, 0, 255])
                        )
                    ),
                    resolution=60,
                    width=groundtrack_width,
                    leadTime=groundtrack_lead_time
                    if groundtrack_lead_time
                    else 100,
                    trailTime=groundtrack_trail_time
                    if groundtrack_trail_time
                    else 100,
                ),
            )
            self.packets.append(pckt)

        self.i += 1

    def add_trajectory(
        self,
        positions,
        epochs,
        groundtrack_show=False,
        groundtrack_lead_time=None,
        groundtrack_trail_time=None,
        groundtrack_width=None,
        groundtrack_color=None,
        id_name=None,
        id_description=None,
        path_width=None,
        path_show=None,
        path_color=None,
        label_fill_color=None,
        label_outline_color=None,
        label_font=None,
        label_text=None,
        label_show=None,
    ):
        """
        Adds trajectory.

        Parameters
        ----------
        positions : ~astropy.coordinates.CartesianRepresentation
            Trajectory to plot.
        epochs : ~astropy.time.Time
            Epochs for positions.
        groundtrack_show : bool
            If set to true, the groundtrack is
            displayed.
        groundtrack_lead_time : float
            The time the animation is ahead of the real-time groundtrack
        groundtrack_trail_time : float
            The time the animation is behind the real-time groundtrack
        groundtrack_width : int
            Groundtrack width
        groundtrack_color : list[int]
            Rgba groundtrack color. By default, it is set to the path color
        id_name : str
            Set orbit name
        id_description : str
            Set orbit description
        path_width : int
            Path width
        path_show : bool
            Indicates whether the path is visible
        path_color : list[int]
            Rgba path color
        label_fill_color : list[int]
            Fill Color in rgba format
        label_outline_color : list[int]
            Outline Color in rgba format
        label_font : str
            Set label font style and size (CSS syntax)
        label_text : str
            Set label text
        label_show : bool
            Indicates whether the label is visible

        """
        positions = (
            positions.represent_as(CartesianRepresentation)
            .get_xyz(1)
            .to_value(u.m)
        )

        epochs = Time(epochs, format="isot")

        if len(epochs) != len(positions):
            raise ValueError("Number of Points and Epochs must be equal.")

        epochs = np.fromiter(
            map(lambda epoch: (epoch - epochs[0]).to_value(u.s), epochs),
            dtype=float,
        )

        positions = np.around(
            np.concatenate([epochs[..., None], positions], axis=1).ravel(), 1
        ).tolist()

        self.trajectories.append([positions, None, label_text, path_color])

        start_epoch = Time(self.start_epoch, format="isot")

        pckt = Packet(
            id=self.i,
            name=id_name,
            description=id_description,
            availability=TimeInterval(
                start=self.start_epoch, end=self.end_epoch
            ),
            position=Position(
                interpolationDegree=5,
                interpolationAlgorithm=InterpolationAlgorithms.LAGRANGE,
                referenceFrame=ReferenceFrames.INERTIAL,
                cartesian=positions,
                # Use explicit UTC timezone, rather than the default, which is a local timezone.
                epoch=start_epoch.datetime.replace(tzinfo=timezone.utc),
            ),
            path=Path(
                show=path_show,
                width=path_width,
                material=Material(
                    solidColor=SolidColorMaterial(
                        color=Color.from_list(path_color)
                    )
                )
                if path_color is not None
                else Material(
                    solidColor=SolidColorMaterial(
                        color=Color.from_list([255, 255, 0])
                    )
                ),
                resolution=120,
            ),
            label=Label(
                text=label_text,
                font=label_font
                if label_font is not None
                else "11pt Lucida Console",
                show=label_show,
                fillColor=Color(rgba=label_fill_color)
                if label_fill_color is not None
                else Color(rgba=[255, 255, 0, 255]),
                outlineColor=Color(rgba=label_outline_color)
                if label_outline_color is not None
                else Color(rgba=[255, 255, 0, 255]),
            ),
            billboard=Billboard(image=PIC_SATELLITE, show=True),
        )

        self.packets.append(pckt)

        if groundtrack_show:
            raise NotImplementedError(
                "Ground tracking for trajectory not implemented yet"
            )

        self.i += 1

    def get_document(self):
        """Retrieves CZML document."""
        return Document(self.packets)
