PIC_SATELLITE = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAADJSURBVDhPnZHRDcMgEEMZjVEYpaNklIzSEfLfD4qNnXAJSFWfhO7w2Zc0Tf9QG2rXrEzSUeZLOGm47WoH95x3Hl3jEgilvDgsOQUTqsNl68ezEwn1vae6lceSEEYvvWNT/Rxc4CXQNGadho1NXoJ+9iaqc2xi2xbt23PJCDIB6TQjOC6Bho/sDy3fBQT8PrVhibU7yBFcEPaRxOoeTwbwByCOYf9VGp1BYI1BA+EeHhmfzKbBoJEQwn1yzUZtyspIQUha85MpkNIXB7GizqDEECsAAAAASUVORK5CYII="

DEFAULTS = {
    "id": "1",
    "name": "NAME",
    "availability": "0000-00-00T00:00:00Z/0000-00-00T00:00:00Z",
    "description": "DESCRIPTION",
    "billboard": {
        "eyeOffset": {"cartesian": [0, 0, 0]},
        "horizontalOrigin": "CENTER",
        "verticalOrigin": "CENTER",
        "image": PIC_SATELLITE,
        "scale": 1,
        "show": True,
        "pixelOffset": {"cartesian2": [0, 0]},
    },
    "label": {
        "fillColor": {"rgba": [255, 255, 0, 255]},
        "font": "11pt Lucida Console",
        "horizontalOrigin": "LEFT",
        "outlineWidth": 2,
        "show": True,
        "style": "FILL_AND_OUTLINE",
        "text": "Text",
        "verticalOrigin": "CENTER",
        "pixelOffset": {"cartesian2": [6, -4]},
        "outlineColor": {"rgba": [2, 100, 0, 255]},
    },
    "path": {
        "show": {
            "interval": "0000-00-00T00:00:00Z/0000-00-00T00:00:00Z",
            "boolean": True,
        },
        "width": 1,
        "resolution": 120,
        "material": {"solidColor": {"color": {"rgba": [255, 255, 0, 255]}}},
    },
}
