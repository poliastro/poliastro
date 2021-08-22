import numpy as np

# Inspired by https://astronomiac.com/color-of-each-planet-in-the-solarsystem/
BODY_COLORS = {
    "Sun": "#ffcc00",
    "Mercury": "#8c8680",
    "Venus": "#e6db67",
    "Earth": "#2a7bd1",
    "Moon": "#999999",
    "Mars": "#cc653f",
    "Jupiter": "#bf8f5c",
    "Saturn": "#decf83",
    "Uranus": "#7ebec2",
    "Neptune": "#3b66d4",
}


def generate_label(epoch, label):
    label_ = f"{epoch.to_value('iso', subfmt='date_hm')}"
    if label:
        label_ += f" ({label})"

    return label_


def generate_sphere(radius, center, num=20):
    u1 = np.linspace(0, 2 * np.pi, num)
    v1 = u1.copy()
    uu, vv = np.meshgrid(u1, v1)

    x_center, y_center, z_center = center

    xx = x_center + radius * np.cos(uu) * np.sin(vv)
    yy = y_center + radius * np.sin(uu) * np.sin(vv)
    zz = z_center + radius * np.cos(vv)

    return xx, yy, zz
