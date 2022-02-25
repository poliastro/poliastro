"""
Created on 10 Nov 2021 22:43:38 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Objective: This file contains functions to compute various characteristic quantities
           for a Cirular Restircted Three Body Problem (CR3BP) model setup

           Features:
               1. Compute CR3BP 'mu', mass ratio of P1-P2 primary bodies
               2. Compute CR3BP 'l*', characterisitic lenght of P1-P2 primary bodies
               3. Comptue CR3BP 't*', charcterisitic time of P1-P2 primary bodies

Useful ideas:
1. To non-dimensioanlize a physical position vector: divide by l*
2. To dimensioanlize a [nd] position vector: multiply by l*
3. To non-dimensionalize a physical veloctiy vecotr, multiply by t*/l*
4. To dimensionalize a [nd] veloctiy vecotr, multiply by l*/t*
5. Position of P1 = [-mu, 0, 0] [nd] {x,y,z}
5. Position of P2 = [1-mu, 0, 0] [nd] {x,y,z}

Physical Qunatities obtained from: JPL’s ephemerides file de405.spk and https://ssd.jpl.nasa.gov/?planet_pos
"""
from poliastro.bodies import (
    Earth,
    Jupiter,
    Mars,
    Mercury,
    Moon,
    Neptune,
    Pluto,
    Saturn,
    Sun,
    Uranus,
    Venus,
)


def mu_calc(mu_pi):
    """Calculate mu of CR3BP
    Parameters
    ----------
    mu_pi : ndarray, float
        mu_pi[0] = mu of P1
        mu_pi[1] = mu of P2

    Returns
    -------
    mu: float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    """
    return mu_pi[1] / (mu_pi[0] + mu_pi[1])


def tstar_calc(mu_pi, dist):
    """Calculate t* of CR3BP
    Parameters
    ----------
    dist: float, [nd]
        Non-dimensional distance between P1 and P2
    mu_pi: ndarray, float
        mu_pi[0] = mu of P1
        mu_pi[1] = mu of P2

    Returns
    -------
    t*: float, [nd]
        sqrt(dist^3/m*)
        Non-dimensional time of P1-P2 system
    """
    return (dist**3 / (mu_pi[0] + mu_pi[1])) ** 0.5


def bodies_char(b1, b2):
    """Calculates mu, dist[km], t* [nd] of the 'body_i - body_j system'
    Parameters
    ----------
    b1 : string
       'body_i'
    b2 : string
       'body_j'

    Returns
    -------
    mu: float
        mu2/(mu1+mu2), mu2<mu1
    dist: float, [nd]
        Non-dimensional distance between P1 and P2
    tstar: float, [nd]
        sqrt(dist^3/m*)
        Non-dimensional time of P1-P2 system
    """
    # Body values
    mu = {}  # km^3 kg^-1 s^-2
    mu["Sun"] = Sun.k.value * 1e-9
    mu["Mercury"] = Mercury.k.value * 1e-9
    mu["Venus"] = Venus.k.value * 1e-9
    mu["Moon"] = Moon.k.value * 1e-9
    mu["Earth"] = Earth.k.value * 1e-9

    mu["Mars"] = Mars.k.value * 1e-9
    mu["Jupiter"] = Jupiter.k.value * 1e-9
    mu["Saturn"] = Saturn.k.value * 1e-9
    mu["Uranus"] = Uranus.k.value * 1e-9
    mu["Neptune"] = Neptune.k.value * 1e-9
    mu["Pluto"] = Pluto.k.value * 1e-9

    mu["Phobos"] = 0.0007112  # Phobos, GM
    mu["Titan"] = 8978.1382  # Titan, GM
    mu["Ganymede"] = 9887.834  # Ganymede, GM
    mu["Titania"] = 228.2  # Titania, GM
    mu["Triton"] = 1427.598  # Triton, GM
    mu["Charon"] = 102.30  # Charon, GM

    #############
    distances = {}  # km
    distances["EarthMoon"] = 384400
    distances["SunEarth"] = 149600000

    distances["SunMars"] = 227944135
    distances["SunJupiter"] = 778279959
    distances["SunSaturn"] = 1427387908
    distances["SunUranus"] = 2870480873
    distances["SunNeptune"] = 4498337290
    distances["SunPluto"] = 5907150229

    distances["MarsPhobos"] = 9376
    distances["JupiterGanymede"] = 1070400
    distances["SaturnTitan"] = 1221865
    distances["UranusTitania"] = 436300
    distances["NeptuneTriton"] = 354759
    distances["PlutoCharon"] = 17536
    #########

    mu_pi = []

    try:
        temp1 = mu[b1]
        temp2 = mu[b2]
        if temp1 == temp2:
            print(
                "Same bodies passed as P1 and P2. Please pass the secondary body of P1 as P2"
            )
            return 0, 0, 0
    except KeyError:
        print("KeyError-> Incorrect/Does Not Exist Input Bodies name")
        return 0, 0, 0

    # Sort P1 and P2 based on their mu values
    if mu[b1] >= mu[b2]:  # b1 is the bigger primary
        mu_pi.append(mu[b1])
        mu_pi.append(mu[b2])
        bodies = b1 + b2
    else:  # b2 is the bigger primary
        mu_pi.append(mu[b2])
        mu_pi.append(mu[b1])

        # To create string to get distance
        bodies = b2 + b1

    mu = mu_calc(mu_pi)

    try:
        dist = distances[bodies]
    except KeyError:
        print(
            "KeyError-> Error in combination of bodies P1-P2, typo/DNE/combo not created"
        )
        return 0, 0, 0
    else:
        tstar = tstar_calc(mu_pi, dist)
        return mu, dist, tstar
