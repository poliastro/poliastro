import pytest  # isort:skip
pytest.importorskip("czml3")

from astropy import units as u

from poliastro.czml.extract_czml import CZMLExtractor
from poliastro.examples import iss, molniya


def test_czml_custom_packet():
    start_epoch = iss.epoch
    end_epoch = iss.epoch + molniya.period

    sample_points = 10

    ellipsoidr = [6373100, 6373100, 6373100]
    pr_map_url = (
        "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg"
    )
    expected_packet = """{
    "id": "custom_properties",
    "properties": {
        "custom_attractor": true,
        "ellipsoid": [
            {
                "array": [
                    6373100,
                    6373100,
                    6373100
                ]
            }
        ],
        "map_url": "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg"
    }
}"""

    extractor = CZMLExtractor(
        start_epoch, end_epoch, sample_points, ellipsoid=ellipsoidr, pr_map=pr_map_url
    )

    pckt = extractor.packets[-1]

    # Test that custom packet parameters where set correctly
    assert repr(pckt) == expected_packet


def test_czml_add_orbit():
    start_epoch = iss.epoch
    end_epoch = iss.epoch + molniya.period

    sample_points = 10

    expected_doc = """[{
    "id": "document",
    "version": "1.0",
    "name": "document_packet",
    "clock": {
        "interval": "2013-03-18T12:00:00.000/2013-03-18T23:59:35.108",
        "currentTime": "2013-03-18T12:00:00.000",
        "multiplier": 60,
        "range": "LOOP_STOP",
        "step": "SYSTEM_CLOCK_MULTIPLIER"
    }
}, {
    "id": "custom_properties",
    "properties": {
        "custom_attractor": true,
        "ellipsoid": [
            {
                "array": [
                    6378137.0,
                    6378137.0,
                    6356752.314245179
                ]
            }
        ],
        "map_url": [
            "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg"
        ]
    }
}, {
    "id": 0,
    "availability": "2013-03-18T10:00:00Z/2013-03-18T21:59:35Z",
    "position": {
        "epoch": "2013-03-18T12:00:00.000",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            -5874061.773508082,
            20159787.8725738,
            40258166.12020527,
            4317.51082821455,
            -11640044.418811519,
            17864068.260633577,
            35673719.98976574,
            8635.0216564291,
            -16129541.019882346,
            13690511.423667355,
            27339319.57262219,
            12952.53248464365,
            -17354822.80901046,
            6974596.3020939445,
            13927946.976716857,
            17270.0433128582,
            -2602428.8843514095,
            -2846586.9793690084,
            -5684502.843750414,
            21587.554141072753,
            17029383.601313725,
            5939060.590166094,
            11860029.944171509,
            25905.0649692873,
            16513858.106876884,
            13042152.675663127,
            26044577.070983447,
            30222.57579750185,
            12264870.944858141,
            17472285.355057023,
            34891347.60591666,
            34540.0866257164,
            6599717.435895659,
            19973508.888030678,
            39886175.58380472,
            38857.59745393095,
            377235.70825277903,
            20840447.559148617,
            41617412.10566203,
            43175.108282145506,
            -5874061.773508167,
            20159787.872573778,
            40258166.12020523,
            47492.61911036005,
            -11640044.418811578,
            17864068.260633536,
            35673719.98976566
        ]
    },
    "billboard": {
        "show": true,
        "image": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVDhPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47WoH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9iaqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYII="
    },
    "label": {
        "text": "Molniya",
        "font": "11pt Lucida Console",
        "style": "FILL",
        "fillColor": {
            "rgba": [
                125,
                80,
                120,
                255
            ]
        },
        "outlineColor": {
            "rgba": [
                255,
                255,
                0,
                255
            ]
        },
        "outlineWidth": 1.0
    },
    "path": {
        "resolution": 120,
        "material": {
            "solidColor": {
                "color": {
                    "rgba": [
                        255,
                        255,
                        0,
                        255
                    ]
                }
            }
        }
    }
}, {
    "id": 1,
    "availability": "2013-03-18T10:00:00Z/2013-03-18T21:59:35Z",
    "position": {
        "epoch": "2013-03-18T12:00:00.000",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            859072.560000001,
            -4137203.6799999997,
            5295568.71,
            4317.51082821455,
            -6275612.3000529,
            -2500051.006459549,
            496823.49381668516,
            8635.0216564291,
            -2951952.719204249,
            3308039.2744516623,
            -5135374.5943018235,
            12952.53248464365,
            5284474.558295589,
            3626657.269700035,
            -2240110.0347697227,
            17270.0433128582,
            4755661.276041964,
            -2067765.0103719088,
            4367882.921366112,
            21587.554141072753,
            -3674993.3280015104,
            -4315653.729708832,
            3705413.3249415695,
            25905.0649692873,
            -5976376.971031915,
            633013.0666635048,
            -3135590.7655361253,
            30222.57579750185,
            1672241.9629143102,
            4536136.094826675,
            -4766186.325441889,
            34540.0866257164,
            6551964.187353159,
            915641.8629497791,
            1510256.5228245915,
            38857.59745393095,
            554555.781288766,
            -4218442.360589035,
            5271693.507946819,
            43175.108282145506,
            -6358691.200073445,
            -2321273.844249726,
            249707.6562815123,
            47492.61911036005,
            -2675736.8087157304,
            3448931.408246476,
            -5194186.303786017
        ]
    },
    "billboard": {
        "show": true,
        "image": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVDhPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47WoH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9iaqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYII="
    },
    "label": {
        "text": "ISS",
        "font": "11pt Lucida Console",
        "style": "FILL",
        "fillColor": {
            "rgba": [
                255,
                255,
                0,
                255
            ]
        },
        "outlineColor": {
            "rgba": [
                255,
                255,
                0,
                255
            ]
        },
        "outlineWidth": 1.0
    },
    "path": {
        "show": false,
        "resolution": 120,
        "material": {
            "solidColor": {
                "color": {
                    "rgba": [
                        255,
                        255,
                        0,
                        255
                    ]
                }
            }
        }
    }
}]"""
    extractor = CZMLExtractor(start_epoch, end_epoch, sample_points)

    extractor.add_orbit(
        molniya, label_text="Molniya", label_fill_color=[125, 80, 120, 255]
    )
    extractor.add_orbit(iss, label_text="ISS", path_show=False)

    assert repr(extractor.packets) == expected_doc


def test_czml_ground_station():
    start_epoch = iss.epoch
    end_epoch = iss.epoch + molniya.period

    sample_points = 10

    expected_doc = """[{
    "id": "document",
    "version": "1.0",
    "name": "document_packet",
    "clock": {
        "interval": "2013-03-18T12:00:00.000/2013-03-18T23:59:35.108",
        "currentTime": "2013-03-18T12:00:00.000",
        "multiplier": 60,
        "range": "LOOP_STOP",
        "step": "SYSTEM_CLOCK_MULTIPLIER"
    }
}, {
    "id": "custom_properties",
    "properties": {
        "custom_attractor": true,
        "ellipsoid": [
            {
                "array": [
                    6378137.0,
                    6378137.0,
                    6356752.314245179
                ]
            }
        ],
        "map_url": [
            "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg"
        ]
    }
}, {
    "id": "GS0",
    "availability": "2013-03-18T10:00:00Z/2013-03-18T21:59:35Z",
    "position": {
        "cartesian": [
            2539356.1623202674,
            4775834.339416022,
            3379897.6662185807
        ]
    },
    "billboard": {
        "show": true,
        "image": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACvSURBVDhPrZDRDcMgDAU9GqN0lIzijw6SUbJJygUeNQgSqepJTyHG91LVVpwDdfxM3T9TSl1EXZvDwii471fivK73cBFFQNTT/d2KoGpfGOpSIkhUpgUMxq9DFEsWv4IXhlyCnhBFnZcFEEuYqbiUlNwWgMTdrZ3JbQFoEVG53rd8ztG9aPJMnBUQf/VFraBJeWnLS0RfjbKyLJA8FkT5seDYS1Qwyv8t0B/5C2ZmH2/eTGNNBgMmAAAAAElFTkSuQmCC"
    },
    "label": {
        "show": true,
        "text": "GS test",
        "font": "11pt Lucida Console",
        "style": "FILL",
        "fillColor": {
            "rgba": [
                120,
                120,
                120,
                255
            ]
        },
        "outlineWidth": 1.0
    }
}, {
    "id": "GS1",
    "availability": "2013-03-18T10:00:00Z/2013-03-18T21:59:35Z",
    "position": {
        "cartesian": [
            4456924.997008477,
            1886774.8000006324,
            4154098.219336245
        ]
    },
    "billboard": {
        "show": true,
        "image": "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACvSURBVDhPrZDRDcMgDAU9GqN0lIzijw6SUbJJygUeNQgSqepJTyHG91LVVpwDdfxM3T9TSl1EXZvDwii471fivK73cBFFQNTT/d2KoGpfGOpSIkhUpgUMxq9DFEsWv4IXhlyCnhBFnZcFEEuYqbiUlNwWgMTdrZ3JbQFoEVG53rd8ztG9aPJMnBUQf/VFraBJeWnLS0RfjbKyLJA8FkT5seDYS1Qwyv8t0B/5C2ZmH2/eTGNNBgMmAAAAAElFTkSuQmCC"
    },
    "label": {
        "show": false,
        "font": "11pt Lucida Console",
        "style": "FILL",
        "outlineWidth": 1.0
    }
}]"""
    extractor = CZMLExtractor(start_epoch, end_epoch, sample_points)

    extractor.add_ground_station(
        [32 * u.degree, 62 * u.degree],
        label_fill_color=[120, 120, 120, 255],
        label_text="GS test",
    )

    extractor.add_ground_station([0.70930 * u.rad, 0.40046 * u.rad], label_show=False)

    assert repr(extractor.packets) == expected_doc


def test_czml_invalid_orbit_epoch_error():
    start_epoch = molniya.epoch
    end_epoch = molniya.epoch + molniya.period

    extractor = CZMLExtractor(start_epoch, end_epoch, 10)

    with pytest.raises(ValueError) as excinfo:
        extractor.add_orbit(iss, label_text="ISS", path_show=False)
    assert (
        "ValueError: The orbit's epoch cannot exceed the constructor's ending epoch"
        in excinfo.exconly()
    )
