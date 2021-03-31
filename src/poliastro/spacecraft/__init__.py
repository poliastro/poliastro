from astropy import units as u


class Spacecraft:
    """
    Class to represent a Spacecraft.
    """

    @u.quantity_input(A=u.km ** 2, C_D=u.one, m=u.kg)
    def __init__(self, A, C_D, m, **metadata):
        """
        Constructor

        Parameters
        ----------
        A: ~astropy.units.Quantity
            Area of the spacecraft (km^2)
        C_D: ~astropy.units.Quantity
            Dimensionless drag coefficient ()
        m: ~astropy.units.Quantity
            Mass of the Spacecraft (kg)
        metadata: Dict[object, dict]
            Optional keyword arguments to Spacecraft

        """

        self._A = A
        self._C_D = C_D
        self._m = m
        self._metadata = metadata

    @property
    def A(self):
        """Returns A, the area of the spacecraft"""
        return self._A

    @property
    def C_D(self):
        """Returns C_D, the drag coefficient"""
        return self._C_D

    @property
    def m(self):
        """Returns m, the mass of the spacecraft"""
        return self._m

    @property
    def ballistic_coefficient(self):
        r"""Calculates the Ballistic coefficient (km^2/kg)

        Returns
        -------
        B: ~astropy.units.Quantity
            Ballistic coefficient (km^2/kg)

        Notes
        -----
        Be aware that you may encounter alternative definitions of the Ballistic
        coefficient, such as the reciprocal:

        .. math::

            \frac{m}{C_D A}

        """

        B = self._C_D * (self._A / self._m)
        return B
