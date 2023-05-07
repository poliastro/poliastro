"""@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University."""


def calculate_mu(mu1, mu2):
    """Calculate mu of CR3BP.

    Parameters
    ----------
    mu1: float, km^3*s^-2
        mu of P1
    mu2: float, km^3*s^-2
        mu of P2

    Returns
    -------
    mu: float, dimensionless
        mu2/(mu1+mu2), mu2<mu1
    """
    return mu2 / (mu1 + mu2)


def calculate_tstar(mu1, mu2, lstar):
    r"""Calculate t* of CR3BP.

    Parameters
    ----------
    mu1: float, km^3*s^-2
        mu of P1
    mu2: float, km^3*s^-2
        mu of P2
    lstar: float, km
        Characterisitc length of P1 - P2 system

    Returns
    -------
    tstar: float, sec
        Characterisitc time of P1-P2 system
        Explicitly converted to u.s as mu1 and mu2 may be of different units

    .. math::

        \sqrt{\frac{l*^3}{M1+M2}}
    """
    return (lstar**3 / (mu1 + mu2)) ** 0.5
