"""
Created on 10 Nov 2021 22:43:38 2022
Updaed on 20 Mar 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Objective: This file contains a class with methods to compute various characteristic quantities
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

class sys_chars:
    """ Computes and stores the properties(mu, l*, t*) of a CR3BP system
        Can be extended to compute properties for BCR4BP model, HR4BP model etc 
    """

    def __init__(self, p1='Earth', p2='Moon'):
        """
        Constructor

        Parameters
        ----------
        p1 : string
           'body_i'. The deafult is Earth
        p2 : string
           'body_j'. The deafult is Moon
        """
        
        mu, lstar, tstar = self.bodies_char_compute(p1,p2)
        self._mu = mu
        self._lstar = lstar
        self._tstar = tstar
        self._p1 = p1
        self._p2 = p2
    
    # All the attributes are made private to make the constant and not be mistakenly changed
    @property
    def mu(self):
        """Mass ration of P1-P2 primary bodies in CR3BP"""
        return self._mu
    @property
    def lstar(self):
        """ Characterisitc Length """
        return self._lstar
    @property
    def tstar(self):
        """ Characterisitc Time """
        return self._tstar
    @property
    def p1(self):
        """ Primary body P1 """
        return self._p1
    @property
    def p2(self):
        """ Primary body P2 """
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
        mu: float, M2/(M1+M2)
            M1 and M2 are mass of Primary Bodies and M2<M1
        """
        return mu_pi[1] / (mu_pi[0] + mu_pi[1])
    
    
    def __tstar_calc(self, mu_pi, dist):
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
    
    
    def bodies_char_compute(self, p1, p2):
        """Calculates mu, dist[km], t* [nd] of the 'body_i - body_j system'
        Parameters
        ----------
        p1 : string
           'body_i'
        p2 : string
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
        
        mu_pi = []
        mu_body, distances = self.__bodies_char()
    
        try:
            temp1 = mu_body[p1]
            temp2 = mu_body[p2]
            if temp1 == temp2:
                print(
                    "Same bodies passed as P1 and P2. Please pass the secondary body of P1 as P2"
                )
                return 0, 0, 0
        except KeyError:
            print("KeyError-> Incorrect/Does Not Exist Input Bodies name")
            return 0, 0, 0
    
        # Sort P1 and P2 based on their mu values
        if mu_body[p1] >= mu_body[p2]:  # p1 is the bigger primary
            mu_pi.append(mu_body[p1])
            mu_pi.append(mu_body[p2])
            bodies = p1 + p2
        else:  # p2 is the bigger primary
            mu_pi.append(mu_body[p2])
            mu_pi.append(mu_body[p1])
    
            # To create string to get distance
            bodies = p2 + p1
    
        mu = self.__mu_calc(mu_pi)
    
        try:
            dist = distances[bodies]
        except KeyError:
            print(
                "KeyError-> Error in combination of bodies P1-P2, typo/DNE/combo not created"
            )
            return 0, 0, 0
        else:
            tstar = self.__tstar_calc(mu_pi, dist)
            return mu, dist, tstar


    def __bodies_char(self):
        """Returns mu value of various celestial bodies and distance between P1-P2 systems
    
        Returns
        -------
        mu_body: dict, float
        distances: dict, float
        """
        # Body values, G*M_body
        mu_body = {}  # km^3 kg^-1 s^-2
        mu_body['Sun'] = 132712440017.99
        mu_body['Moon'] = 4902.8005821478
        mu_body['Earth'] = 398600.4415
    
        mu_body['Mars'] = 42828.314258067 # Mars, GM
        mu_body['Jupiter'] = 126712767.8578 # Jupiter, GM
        mu_body['Saturn'] = 37940626.061137 # Saturn, GM
        mu_body['Uranus'] = 5794549.0070719 # Uranus, GM
        mu_body['Neptune'] = 6836534.0638793 # Neptune, GM
        mu_body['Pluto'] = 981.600887707 # Pluto, GM
        
        mu_body["Phobos"] = 0.0007112  # Phobos, GM
        mu_body["Titan"] = 8978.1382  # Titan, GM
        mu_body["Ganymede"] = 9887.834  # Ganymede, GM
        mu_body["Titania"] = 228.2  # Titania, GM
        mu_body["Triton"] = 1427.598  # Triton, GM
        mu_body["Charon"] = 102.30  # Charon, GM
    
        #########
        distances = {}  # km, diistance between the two primaries
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
        
        return mu_body, distances