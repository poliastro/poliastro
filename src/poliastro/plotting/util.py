import numpy as np

BODY_COLORS = {"Sun": "#ffcc00", "Earth": "#204a87", "Jupiter": "#ba3821"}


def generate_label(orbit, label):
    epoch = orbit.epoch.copy()
    epoch.out_subfmt = "date_hm"
    label_ = f"{epoch.iso}"
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


def generate_circle(radius, center, num=500):
    u1 = np.linspace(0, 2 * np.pi, num)
    x_center, y_center, z_center = center

    xx = x_center + radius * np.cos(u1)
    yy = y_center + radius * np.sin(u1)

    return xx, yy
