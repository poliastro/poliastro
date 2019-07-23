import pytest  # noqa: E402 isort:skip

from astropy import units as u  # noqa: E402

from poliastro.czml.extract_czml import CZMLExtractor  # noqa: E402
from poliastro.examples import iss, molniya  # noqa: E402


def test_czml_custom_packet():
    start_epoch = iss.epoch
    end_epoch = iss.epoch + molniya.period

    sample_points = 10

    ellipsoidr = [6373100, 6373100, 6373100]
    pr_map_url = (
        "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg"
    )
    scene = False
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
        "map_url": "https://upload.wikimedia.org/wikipedia/commons/c/c4/Earthmap1000x500compac.jpg",
        "scene3D": false
    }
}"""

    extractor = CZMLExtractor(
        start_epoch,
        end_epoch,
        sample_points,
        ellipsoid=ellipsoidr,
        pr_map=pr_map_url,
        scene3D=scene,
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
        ],
        "scene3D": true
    }
}, {
    "id": 0,
    "availability": "2013-03-18T12:00:00Z/2013-03-18T23:59:35Z",
    "position": {
        "epoch": "2013-03-18T12:00:00.000",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            -5874061.7735,
            20159787.8726,
            40258166.1202,
            4317.5108,
            -11640044.4188,
            17864068.2606,
            35673719.9898,
            8635.0217,
            -16129541.0199,
            13690511.4237,
            27339319.5726,
            12952.5325,
            -17354822.809,
            6974596.3021,
            13927946.9767,
            17270.0433,
            -2602428.8844,
            -2846586.9794,
            -5684502.8438,
            21587.5541,
            17029383.6013,
            5939060.5902,
            11860029.9442,
            25905.065,
            16513858.1069,
            13042152.6757,
            26044577.071,
            30222.5758,
            12264870.9449,
            17472285.3551,
            34891347.6059,
            34540.0866,
            6599717.4359,
            19973508.888,
            39886175.5838,
            38857.5975,
            377235.7083,
            20840447.5591,
            41617412.1057,
            43175.1083,
            -5874061.7735,
            20159787.8726,
            40258166.1202,
            47492.6191,
            -11640044.4188,
            17864068.2606,
            35673719.9898
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
    "availability": "2013-03-18T12:00:00Z/2013-03-18T23:59:35Z",
    "position": {
        "epoch": "2013-03-18T12:00:00.000",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            859072.56,
            -4137203.68,
            5295568.71,
            4317.5108,
            -6275612.3001,
            -2500051.0065,
            496823.4938,
            8635.0217,
            -2951952.7192,
            3308039.2745,
            -5135374.5943,
            12952.5325,
            5284474.5583,
            3626657.2697,
            -2240110.0348,
            17270.0433,
            4755661.276,
            -2067765.0104,
            4367882.9214,
            21587.5541,
            -3674993.328,
            -4315653.7297,
            3705413.3249,
            25905.065,
            -5976376.971,
            633013.0667,
            -3135590.7655,
            30222.5758,
            1672241.9629,
            4536136.0948,
            -4766186.3254,
            34540.0866,
            6551964.1874,
            915641.8629,
            1510256.5228,
            38857.5975,
            554555.7813,
            -4218442.3606,
            5271693.5079,
            43175.1083,
            -6358691.2001,
            -2321273.8442,
            249707.6563,
            47492.6191,
            -2675736.8087,
            3448931.4082,
            -5194186.3038
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
        molniya, rtol=1e-4, label_text="Molniya", label_fill_color=[125, 80, 120, 255]
    )
    extractor.add_orbit(iss, rtol=1e-4, label_text="ISS", path_show=False)

    assert repr(extractor.packets) == expected_doc


def test_czml_groundtrack():

    start_epoch = molniya.epoch
    end_epoch = molniya.epoch + molniya.period

    sample_points = 10

    expected_doc = """[{
    "id": "document",
    "version": "1.0",
    "name": "document_packet",
    "clock": {
        "interval": "2000-01-01T12:00:00.000/2000-01-01T23:59:35.108",
        "currentTime": "2000-01-01T12:00:00.000",
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
        ],
        "scene3D": true
    }
}, {
    "id": 0,
    "availability": "2000-01-01T12:00:00Z/2000-01-01T23:59:35Z",
    "position": {
        "epoch": "2000-01-01T12:00:00.000",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            10140093.6393,
            -800580.7546,
            -1598722.8246,
            4317.5108,
            17593269.6479,
            8807305.9945,
            17587783.6632,
            8635.0217,
            15276972.8693,
            14840949.3322,
            29636690.9896,
            12952.5325,
            10390291.9946,
            18544487.5609,
            37032485.9349,
            17270.0433,
            4459288.7909,
            20454288.4376,
            40846270.159,
            21587.5541,
            -1827745.3661,
            20778744.495,
            41494193.9337,
            25905.065,
            -7974181.1167,
            19549661.224,
            39039771.3569,
            30222.5758,
            -13412369.429,
            16641263.6853,
            33231835.6785,
            34540.0866,
            -17116767.3657,
            11688710.5829,
            23341815.6655,
            38857.5975,
            -15852598.671,
            3766095.5927,
            7520719.1144,
            43175.1083,
            10140093.6393,
            -800580.7546,
            -1598722.8246,
            47492.6191,
            17593269.6479,
            8807305.9945,
            17587783.6632
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
    "id": "groundtrack0",
    "availability": "2000-01-01T12:00:00Z/2000-01-01T23:59:35Z",
    "position": {
        "epoch": "2000-01-01T12:00:00.000",
        "interpolationAlgorithm": "LAGRANGE",
        "interpolationDegree": 5,
        "referenceFrame": "INERTIAL",
        "cartesian": [
            0.0,
            6280728.650052842,
            -495876.0310202439,
            -990241.6517299684,
            4317.5108,
            4245757.090389717,
            2125453.852544122,
            4244433.16605073,
            8635.0217,
            2663916.0370143135,
            2587884.609889886,
            5167886.049854508,
            12952.5325,
            1548100.6194718895,
            2763034.174623665,
            5517651.648479046,
            17270.0433,
            618022.3743442999,
            2834803.2385849194,
            5660971.157461388,
            21587.5541,
            -250341.4520219223,
            2846010.9516345025,
            5683352.443234859,
            25905.065,
            -1142958.8856769665,
            2802101.1102361814,
            5595666.516165981,
            30222.5758,
            -2159933.1837527766,
            2679915.7933351607,
            5351668.076502109,
            34540.0866,
            -3490752.5841853595,
            2383767.83897665,
            4760274.277201323,
            38857.5975,
            -5630877.25034594,
            1337725.4037362342,
            2671375.803029601,
            43175.1083,
            6280728.650052842,
            -495876.0310202439,
            -990241.6517299684,
            47492.6191,
            4245757.090389717,
            2125453.852544122,
            4244433.16605073
        ]
    },
    "path": {
        "show": true,
        "width": 2,
        "resolution": 60,
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
        molniya,
        rtol=1e-4,
        label_text="Molniya",
        show_groundtrack=True,
        label_fill_color=[125, 80, 120, 255],
    )

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
        ],
        "scene3D": true
    }
}, {
    "id": "GS0",
    "availability": "2013-03-18T12:00:00Z/2013-03-18T23:59:35Z",
    "position": {
        "cartesian": [
            2536956.641087247,
            4771321.49632892,
            3399461.1745119346
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
    "availability": "2013-03-18T12:00:00Z/2013-03-18T23:59:35Z",
    "position": {
        "cartesian": [
            4450567.756252572,
            1884083.5539814981,
            4176129.520641576
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
