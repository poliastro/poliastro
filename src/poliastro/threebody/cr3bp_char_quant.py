"""
@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University

Useful ideas:
1. To non-dimensioanlize a physical position vector: divide by l*
2. To dimensioanlize a [nd] position vector: multiply by l*
3. To non-dimensionalize a physical veloctiy vecotr, multiply by t*/l*
4. To dimensionalize a [nd] veloctiy vecotr, multiply by l*/t*
5. Position of P1 = [-mu, 0, 0] [nd] {x,y,z}
5. Position of P2 = [1-mu, 0, 0] [nd] {x,y,z}
"""
from astropy import units as u

import poliastro.constants.general as general
import poliastro.constants.mean_elements as dist


class SystemChars:
    """Computes and stores the properties (mu, l*, t*) of a CR3BP system
    'mu': mass ratio of P1-P2 primary bodies
    'l*': characterisitic lenght of P1-P2 primary bodies
    't*': characterisitic time of P1-P2 primary bodies
    If P2 is more massive than P1 then swap the Pi, so that P1 is the more massive body
    """

    def __init__(self, p1="Earth", p2="Moon"):
        """
        Constructor

        Parameters
        ----------
        p1 : string
           'body_i'. The deafult is Earth
        p2 : string
           'body_j'. The deafult is Moon
        """

        mu, lstar, tstar, p1, p2 = self.BodiesCharCompute(
            p1.lower(), p2.lower()
        )  # lower() to ensure lowercase characters
        units_nd = u.def_unit(
            "[nd]"
        )  # Define custome units, [nd]: non-dimensional unit

        self._mu = mu * units_nd
        self._lstar = lstar
        self._tstar = tstar
        self._p1 = p1
        self._p2 = p2

    # All the attributes are made private to make them constant and avoid being mistakenly changed
    @property
    def mu(self):
        """Mass ration of P1-P2 primary bodies in CR3BP"""
        return self._mu

    @property
    def lstar(self):
        """Characterisitc Length"""
        return self._lstar

    @property
    def tstar(self):
        """Characterisitc Time"""
        return self._tstar

    @property
    def p1(self):
        """Primary body P1"""
        return self._p1

    @property
    def p2(self):
        """Primary body P2"""
        return self._p2

    def __mu_calc(self, mu_pi):
        """Calculate mu of CR3BP
        Parameters
        ----------
        mu_pi : ndarray, float
            mu_pi[0] = mu of P1
            mu_pi[1] = mu of P2

        Returns
        -------
        mu: float, M2/(M1+M2) = mu2/(mu1+mu2)
            M1 and M2 are mass of Primary Bodies and M2<M1
        """
        return mu_pi[1] / (mu_pi[0] + mu_pi[1])

    def __tstar_calc(self, mu_pi, lstar):
        """Calculate t* of CR3BP
        Parameters
        ----------
        mu_pi: ndarray, float
            mu_pi[0] = mu of P1
            mu_pi[1] = mu of P2
        lstar: float, km
            Characterisitc distance between P1 and P2

        Returns
        -------
        t*: float, sec
            Characterisitc time of P1-P2 system

        .. math::

            \sqrt{\frac{l*^3}{M1+M2}}
        """
        return (lstar**3 / (mu_pi[0] + mu_pi[1])) ** 0.5

    def BodiesCharCompute(self, p1, p2):
        """Calculates mu [nd], dist[km], t* [s] of the 'body_i - body_j system'
        Also, if M2>M1, then swaps p1 and p2, so that M1>M2

        Parameters
        ----------
        p1 : string
           'body_i'
        p2 : string
           'body_j'

        Returns
        -------
        mu: float, [nd]
            M2/(M1+M2) = mu2/(mu1+mu2), mu2<mu1
        lstar: float, km
            Characterisitc length between P1 and P2
        tstar: float, sec
            Characterisitc time of P1-P2 system
        p1 : string
           'body_i', more massive body
        p2 : string
           'body_j', less massive body
        """

        mu_pi = []
        mu_body, distances = self.__bodies_char()

        try:
            temp1 = mu_body[p1]
            temp2 = mu_body[p2]
            if p1 == p2:
                print(
                    "Same bodies passed as P1 and P2. Please pass the secondary body of P1 as P2"
                )
                return 0, 0, 0, temp1, temp2
        except KeyError:
            print(
                "KeyError-> Incorrect/Bodies Does Not Exist Input Bodies name"
            )
            return 0, 0, 0, temp1, temp2

        # Sort P1 and P2 based on their mu values
        if mu_body[p1] >= mu_body[p2]:  # if p1 is the bigger primary
            mu_pi.append(mu_body[p1])
            mu_pi.append(mu_body[p2])
            bodies = p1 + p2
        else:  # if p2 is the bigger primary
            mu_pi.append(mu_body[p2])
            mu_pi.append(mu_body[p1])

            # Create a string to get distance of the P1-P2 system
            bodies = p2 + p1
            # swap p1 and p2, as p1 should be the more massive body
            temp_pi = p2
            p2 = p1
            p1 = temp_pi

        mu = self.__mu_calc(mu_pi)

        try:
            lstar = distances[bodies]
        except KeyError:
            print(
                "KeyError-> Error in combination of bodies P1-P2, typo/DNE/system not created"
            )
            return 0, 0, 0, p1, p2
        else:
            tstar = self.__tstar_calc(mu_pi, lstar)
            return mu, lstar, tstar, p1, p2

    def __bodies_char(self):
        """Returns mass ratio and mean distance of various P1-P2 systems

        Returns
        -------
        mu_body: dict, float, km^3*s^-2
        distances: dict, float, km
        """
        # Body values, G*M_body
        mu_body = {}  # km^3 s^-2
        mu_body["sun"] = general.GM_sun.to(u.km**3 * u.s**-2)
        mu_body["mercury"] = general.GM_mercury.to(u.km**3 * u.s**-2)
        mu_body["venus"] = general.GM_venus.to(u.km**3 * u.s**-2)
        mu_body["earth"] = general.GM_earth.to(u.km**3 * u.s**-2)
        mu_body["mars"] = general.GM_mars.to(u.km**3 * u.s**-2)
        mu_body["jupiter"] = general.GM_jupiter.to(u.km**3 * u.s**-2)
        mu_body["saturn"] = general.GM_saturn.to(u.km**3 * u.s**-2)
        mu_body["uranus"] = general.GM_uranus.to(u.km**3 * u.s**-2)
        mu_body["neptune"] = general.GM_neptune.to(u.km**3 * u.s**-2)
        mu_body["pluto"] = general.GM_pluto.to(u.km**3 * u.s**-2)

        mu_body["moon"] = general.GM_moon.to(u.km**3 * u.s**-2)
        mu_body["phobos"] = general.GM_phobos.to(u.km**3 * u.s**-2)
        mu_body["deimos"] = general.GM_deimos.to(u.km**3 * u.s**-2)
        mu_body["europa"] = general.GM_europa.to(u.km**3 * u.s**-2)
        mu_body["ganymede"] = general.GM_ganymede.to(u.km**3 * u.s**-2)
        mu_body["enceladus"] = general.GM_enceladus.to(u.km**3 * u.s**-2)
        mu_body["titania"] = general.GM_titan.to(u.km**3 * u.s**-2)
        mu_body["triton"] = general.GM_triton.to(u.km**3 * u.s**-2)
        mu_body["charon"] = general.GM_charon.to(u.km**3 * u.s**-2)

        distances = {}  # km, diistance between the two primaries
        distances["sunmercury"] = dist.mean_a_mercury
        distances["sunvenus"] = dist.mean_a_venus
        distances["sunearth"] = dist.mean_a_earth
        distances["sunmars"] = dist.mean_a_mars
        distances["sunjupiter"] = dist.mean_a_jupiter
        distances["sunsaturn"] = dist.mean_a_saturn
        distances["sunuranus"] = dist.mean_a_uranus
        distances["sunneptune"] = dist.mean_a_neptune

        distances["earthmoon"] = dist.mean_a_moon
        distances["marsphobos"] = dist.mean_a_phobos
        distances["marsdeimos"] = dist.mean_a_deimos
        distances["jupitereuropa"] = dist.mean_a_europa
        distances["jupiterganymede"] = dist.mean_a_ganymede
        distances["saturnenceladus"] = dist.mean_a_enceladus
        distances["saturntitan"] = dist.mean_a_titan
        distances["uranustitania"] = dist.mean_a_titania
        distances["neptunetriton"] = dist.mean_a_triton
        distances["plutocharon"] = dist.mean_a_charon

        return mu_body, distances
