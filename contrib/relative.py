# -*- coding: utf-8 -*-

##############################################################################
##############################################################################
##                                                                          ##
##   FILE SUMMARY: RELATIVE ORBIT CLASS FOR PLOTTING AND VISUALISATION      ##
##                                                                          ##
##                                                                          ##
##                                                                          ##
##   FILE DESCRIPTION:                                                      ##
##                                                                          ##
##   This file is a prototype file for a newly proposed "Relative Orbit"    ##
##   object in poliastro. The class 'relative' is defined by 02x instances  ##
##   of the 'twobody' object - one chief, and one deputy spacecraft, and    ##
##   it takes in the classical elements from each of them. The class then   ##
##   automatically computes the states of relative motion by linearizing    ##
##   the equations of motion based on the Clohessy-Wiltshire equations.     ##
##                                                                          ##
##   The user may then either sample the state vectors of the relative      ##
##   motion over time by providing an epoch input, and plot the relative    ##
##   trajectory in a 3D and interactive matplotlib plot.                    ##
##                                                                          ##
##   The accuracy of this algorithm has been validated in AGI's STK 10.     ##
##                                                                          ##
##                                                                          ##
##                                                                          ##
##   BY: SAMUEL Y.W. LOW (sammmlow@gmail.com)                               ##
##                                                                          ##
##   LAST MODIFIED: 02-05-2021 00:53 (GMT)                                  ##
##                                                                          ##
##   IN CONTRIBUTION TO: POLIASTRO 0.15.dev0                                ##
##                                                                          ##
##                                                                          ##
##                                                                          ##
##   REFERENCES:                                                            ##
##                                                                          ##
##   [1] Clohessy, W. H., Wiltshire, R. S. (1960). Terminal guidance        ##
##       system for satellite rendezvous. Journal of the Aerospace          ##
##       Sciences, 27(9), 653-658. doi:10.2514/8.8704                       ##
##                                                                          ##
##   [2] Montenbruck, O., Gill, E. (2013). Satellite Orbits Models,         ##
##       Methods and Applications. Springer, Berlin.                        ##
##                                                                          ##
##   [3] D'Amico, S., &; Montenbruck, O. (2006). Proximity operations of    ##
##       formation-flying spacecraft using an eccentricity/inclination      ##
##       vector separation. Journal of Guidance, Control, and Dynamics,     ##
##       29(3), 554-563. doi:10.2514/1.15114                                ##
##                                                                          ##
##   [4] D'Amico, S. (2010) Autonomous Formation Flying in Low Earth Orbit  ##
##       PhD Dissertation, TU Delft.                                        ##
##                                                                          ##
##                                                                          ##
##                                                                          ##
##############################################################################
##############################################################################

# Import numpy.
import numpy as np

# Import matplotlib (possibly change to plotly for poliastro?)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import the relevant AstroPy libraries.
from astropy import time, units as u

# In order to plot astropy quantities, this library needs to be turned on.
from astropy.visualization import quantity_support

# Import the relevant PoliAstro libraries
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.angles import nu_to_E, E_to_M

# Below is the definition of the relative orbits class (RelativeOrb).
# For the actual example script of how to use the class, scroll to the bottom!

##############################################################################
##############################################################################
##                                                                          ##
##   THE RELATIVE ORBIT CLASS:                                              ##
##   SCROLLDOWN TO THE BOTTOM FOR THE ACTUAL SCRIPT!                        ##
##                                                                          ##
##############################################################################
##############################################################################

# Let us define a relative orbits class that can be defined by two orbits!
class RelativeOrb:
    
    def __init__(self, satC, satD):
        
        """
        
        When this class is initialized, it requires two poliastro twobody
        objects: a chief orbiter (satC) and a deputy orbiter (satD). After
        initialialization, a number of computations must be made.
        
        First, the relative orbit is defined by the separation of the relative
        inclination vector [ix,iy], and relative eccentricity vector [ex,ey].
        
        Second, it is also defined by the difference between the semi-major
        axis. These differential terms that are not time-dependent can be
        initialized at the beginning.
        
        The above are time-independent and thus can be computed before any
        propagation is performed. The time-dependent components of relative
        states (which is actually just 'du') will be initialized but not
        computed until propagation.
        
        Since the equations of motion are realized by linearizing about the
        chief orbiter motion, all relative elements are defined in the order
        of deputy-minus-chief, and absolute orbit elements are taken as the
        chief elements. See D'Amico's paper [3] for more details, although I
        have slightly modified his equations to take into account effects at
        low-inclination orbits that was not captured in his paper.
        
        Computed parameters
        -------------------
        - ix --> Relative inclination vector x-component
        - iy --> Relative inclination vector y-component
        - ex --> Relative eccentricity vector x-component
        - ey --> Relative eccentricity vector y-component
        - da --> Relative semi-major axis (normalized over the Chief SMA)
        - dR --> Relative in-track separation due to RAAN differential
        - du --> Relative argument of latitude (time dependent, default zero)
        
        """
        
        # Initialize the chief and deputy satellites that were class inputs.
        self.satC = satC
        self.satD = satD
        
        # Below are all the auxiliary parameters used in the computation of
        # the state transition matrix that converts the time-dependent anomaly
        # values into relative position and velocity values in the VVLH frame.
        # Ensure that all the values below are of the same type (numpy.float64)
        self.ix = (satD.inc - satC.inc).to_value(u.rad)
        self.iy = (np.sin(satC.inc) * (satD.raan - satC.raan)).to_value(u.rad)
        self.ex = (satD.ecc*np.cos(satD.argp)) - (satC.ecc*np.cos(satC.argp))
        self.ex = self.ex.to_value(u.one)
        self.ey = (satD.ecc*np.sin(satD.argp)) - (satC.ecc*np.sin(satC.argp))
        self.ey = self.ey.to_value(u.one)
        self.da = ((satD.a - satC.a) /  satC.a).to_value(u.one)
        self.dR = ((satD.raan - satC.raan) * np.cos(satC.inc)).to_value(u.rad)
        
        # With the above, we can already define the state transition matrix
        # from classical orbit elements to VVLH frame relative state vectors.
        self.M = np.zeros((6,4))
        
        # Column 1 of Linearized Relative EOM Matrix
        self.M[0][0] = self.da
        self.M[1][0] = self.dR # This element will be updated in the for loop.
        self.M[2][0] = 0.0
        self.M[3][0] = 0.0
        self.M[4][0] = self.da * (-1.5)
        self.M[5][0] = 0.0
        
        # Column 2 of Linearized Relative EOM Matrix
        self.M[0][1] = 0.0
        self.M[1][1] = self.da * (-1.5)
        self.M[2][1] = 0.0
        self.M[3][1] = 0.0
        self.M[4][1] = 0.0
        self.M[5][1] = 0.0
        
        # Column 3 of Linearized Relative EOM Matrix
        self.M[0][2] = self.ex * (-1)
        self.M[1][2] = 0.0   # This was mistakenly = -2*ey in the paper
        self.M[2][2] = self.iy * (-1)
        self.M[3][2] = self.ey * (-1)
        self.M[4][2] = 0.0   # This was mistakenly = 2*ex in the paper
        self.M[5][2] = self.ix
        
        # Column 4 of Linearized Relative EOM Matrix
        self.M[0][3] = self.ey * (-1)
        self.M[1][3] = 0.0   # This was mistakenly = 2*ex in the paper
        self.M[2][3] = self.ix
        self.M[3][3] = self.ex
        self.M[4][3] = 0.0   # This was mistakenly = 2*ey in the paper
        self.M[5][3] = self.iy
        
        # Finally, define the output of the relative trajectory propagation.
        self.relPosArray = np.array([]).astype( type(1*u.km) )
        self.relVelArray = np.array([]).astype( type(1*u.km/u.s) )
    
    # Method to return a list of the computed eccentricity vector separation.
    def get_eccentricity_separation(self):
        
        """
        Returns the eccentricity separation vector [ex, ey]
        """
        
        return [ self.ex, self.ey ]
    
    # Method to return a list of the computed inclination vector separation.
    def get_inclination_separation(self):
        
        """
        Returns the inclination separation vector [ex, ey]
        """
        
        return [ self.ix, self.iy ]
    
    # Method to return direction cosine matrix for a rotation t about X-axis
    def _dcmX(self, t):
        
        """
        Input theta is the scalar angle (IN RADIANS).
        Output is a 3x3 direction cosine matrix.
        """
        
        dcm = np.array([[ 1.0,    0.0,       0.0       ],
                        [ 0.0,    np.cos(t), np.sin(t) ],
                        [ 0.0, -1*np.sin(t), np.cos(t) ]])
    
        return dcm
    
    # Method to return direction cosine matrix for a rotation t about Z-axis
    def _dcmZ(self, t):
    
        """
        Input theta is the scalar angle (IN RADIANS).
        Output is a 3x3 direction cosine matrix.
        """
        
        dcm = np.array([[    np.cos(t), np.sin(t), 0.0 ],
                        [ -1*np.sin(t), np.cos(t), 0.0 ],
                        [    0.0,       0.0,       1.0 ]])
        
        return dcm
    
    # Method to solve for Kepler's equation with mean anomaly input
    def _solve_kepler(self, M, ecc):
        
        """
        Input a float mean anomaly (rad) and eccentricity, and solves for
        the eccentric anomaly.
        """
        
        E1 = M.to_value(u.rad) # Initialise eccentric anomaly
        e = ecc.to_value(u.one) # Initialise the float eccentricity
        residual = 1.0 # Initialise convergence residual
        
        while residual >= 0.000001:
            
            fn = E1 - (e*np.sin(E1)) - M.to_value(u.rad)
            fd = 1 - (e*np.cos(E1))
            E2 = E1 - (fn/fd)
            residual = abs(E2-E1) # Compute residual
            E1 = E2 # Update the eccentric anomaly
            
        return E2 * u.rad
    
    # Method to solve for the orbit position, velocity and true anomaly.
    # (Requires the _dcmZ, _dcmX, and _solve_kepler internal methods)
    def _solve_posn(self, a, e, i, w, R, M, mu):
        
        """
        Inputs: Keplerian elements and gravitational constant (Astropy units).
                - a    -> Semi-major axis (km)
                - e    -> Eccentricity (unit-less)
                - i    -> Inclination (degrees)
                - w    -> Argument of Perigee (degrees)
                - R    -> Right Angle of Asc Node (degrees)
                - M    -> Mean Anomaly (degrees)
                - mu   -> Gravitational constant (km^3/s^2)
        
        Output: A 1x3 position and 1x3 velocity vector in an inertial frame.
        """
        
        mu = mu.to(u.km**3/u.s**2)
        
        # The general flow of the program, is to first solve for the radial
        # position and velocity (in the inertial frame) via Kepler's equation.
        # Thereafter, we obtain the inertial coordinates in the Hill frame,
        # by performing a 3-1-3 Euler Angle rotation using an appropriate DCM.
        
        # First, let us solve for the eccentric anomaly.
        eccAnom = self._solve_kepler(M,e)
        
        # With the eccentric anomaly, we can solve for position and velocity
        # in the local orbital frame, using the polar equation for an ellipse.
        pos_X = a * ( np.cos(eccAnom) - e)
        pos_Y = a * np.sqrt( 1 - e**2 ) * np.sin(eccAnom)
        pos_norm = np.sqrt( pos_X**2 + pos_Y**2 )
        vel_const = np.sqrt( mu * a ) / pos_norm
        vel_X = vel_const * ( -1 * np.sin(eccAnom) )
        vel_Y = vel_const * ( np.sqrt( 1 - e**2 ) * np.cos(eccAnom) )
        
        # To perform the conversion from local orbit plane to an ECI frame, we
        # need perform the 313 Euler angle rotation in the following sequence:
        # Right Angle of Ascending Node -> Inclination -> Argument of Latitude.
        # Now, let us get us the DCM that converts to the hill-frame.
        DCM_HN = np.matmul( self._dcmZ(w),
                            np.matmul( self._dcmX(i), 
                                       self._dcmZ(R) ) )
        
        # Notice that the hill frame computation does not include a rotation
        # of the true anomaly, and that's because the true anomaly has already
        # been accounted for when computing pos_X and pos_Y using information 
        # from the eccentric anomaly. Including true anomaly in the DCM 
        # rotation would double-count that anomaly rotation.
        
        # The current coordinates are in the local hill frame, and thus 
        # conversion from hill to inertial would be the transpose of HN.
        DCM_NH = np.transpose(DCM_HN)
        
        # For the matrix multiplication below, we will 
        
        # With the hill frame, we can now convert it to the ECI frame.
        pos = np.matmul(DCM_NH, np.array([ pos_X.to_value(u.km),
                                           pos_Y.to_value(u.km),
                                           0.0 ]))
        vel = np.matmul(DCM_NH, np.array([ vel_X.to_value(u.km/u.s),
                                           vel_Y.to_value(u.km/u.s),
                                           0.0 ]))
        
        # Finally, let us not forget to compute the true anomaly.
        trueAnom = np.arctan2( pos_Y.to_value(u.km), pos_X.to_value(u.km) ) 
        
        # Position vector 1x3 (m), velocity vetor 1x3 (m), true anomaly (rad)
        return pos * u.km, vel * (u.km/u.s), trueAnom * u.rad
        
    # Propagate method, that must be called in order to store values for plots.
    def propagate(self, duration=43200, step=60):
        
        """
        This is the core method used for the relative trajectory generation.
        Propagates the relative orbit trajectory for one full orbit period.

        Input:
        - duration - integer number of seconds (optional, default = 12 hours)
        - step - integer time step size (optional, default = 1 minute)
        
        Returns two numpy arrays
        - 1x3 position vector
        - 1x3 velocity vector
        - float true anomaly (rad)
        """
        
        # Check if the orbit is an ellipse (closed)
        if self.satC.ecc < 1 and self.satD.ecc < 1:
            
            # Get the initial mean anomaly of the chief.
            mC = E_to_M( nu_to_E( self.satC.nu, self.satC.ecc ), self.satC.ecc)
            
            # Get the mean anomaly of the deputy.
            mD = E_to_M( nu_to_E( self.satD.nu, self.satD.ecc ), self.satD.ecc)
            
            # Get the mean motion of the chief.
            nC = self.satC.n
            
            # Get the mean motion of the deputy.
            nD = self.satD.n
            
            # Initialise relative position and velocity component arrays.
            relPosArrayX, relPosArrayY, relPosArrayZ = [], [], []
            relVelArrayX, relVelArrayY, relVelArrayZ = [], [], []
            
            # Get the gravitational constant in u.km**3 / u.s**2
            mu = self.satC.attractor.k.to( u.km**3 / u.s**2 )
            
            # Initialise pi in terms of astropy units
            pi = np.pi * 1 * u.rad
            
            # For each sample...
            for t in range( 0, duration, step ):
                
                # Get a time step in Astropy units (seconds)
                ts = step * u.s
                
                # Update the mean anomaly of the chief.
                mC = ( ( mC + pi + ( nC * ts ) ) % ( 2 * pi ) ) - pi
                
                # Update the mean anomaly of the deputy.
                mD = ( ( mD + pi + ( nD * ts ) ) % ( 2 * pi ) ) - pi
                
                # Compute the chief position, velocity and true anomaly.
                pC, vC, nuC = self._solve_posn(self.satC.a, 
                                               self.satC.ecc, 
                                               self.satC.inc, 
                                               self.satC.argp, 
                                               self.satC.raan, 
                                               mC, mu)
                
                # Compute the deputy position, velocity and true anomaly.
                pD, vD, nuD = self._solve_posn(self.satD.a, 
                                               self.satD.ecc, 
                                               self.satD.inc, 
                                               self.satD.argp, 
                                               self.satD.raan, 
                                               mD, mu)
                
                # Get the argument of latitude of the chief.
                uC = nuC + self.satC.argp
                uC = ( uC + pi ) % ( 2 * pi ) - pi
                
                # Get the argument of latitude of the deputy.
                uD = nuD + self.satD.argp
                uD = ( uD + pi ) % ( 2 * pi ) - pi
                
                # Get the relative argument of latitude.
                du = uD - uC
                du = ( du + pi ) % ( 2 * pi ) - pi
                
                # Save the chief initial argument of latitude.
                if t == 0:
                    uC0 = uC
                
                # Compute the deputy elapsed argument of latitude.
                uD_elapsed = uD - uC0
                uD_elapsed = ( uD_elapsed + pi ) % ( 2 * pi ) - pi
                uD_elapsed = uD_elapsed.to_value(u.rad)
                
                # Compute the velocity magnitude
                vCMag = np.sqrt( vC[0]**2 + vC[1]**2 + vC[2]**2 )
                
                # Initialize the time-dependent input vector
                uVect = np.array([ 1.0, uD_elapsed, np.cos(uC), np.sin(uC) ])
                
                # Update Row 1 Column 0 of the state transition matrix.
                self.M[1][0] = du.to_value(u.rad) + self.dR
                
                # We may now compute the normalized relative state vectors
                relPos = np.matmul( self.M[:3], uVect )
                relVel = np.matmul( self.M[3:], uVect )
                
                # Un-normalize the relative position vectors
                relPosArrayX.append( (relPos[0] * self.satC.a).to_value(u.km) )
                relPosArrayY.append( (relPos[1] * self.satC.a).to_value(u.km) )
                relPosArrayZ.append( (relPos[2] * self.satC.a).to_value(u.km) )
                
                # Un-normalize the relative velocity vectors
                relVelArrayX.append( (relVel[0] * vCMag).to_value(u.km/u.s) )
                relVelArrayY.append( (relVel[1] * vCMag).to_value(u.km/u.s) )
                relVelArrayZ.append( (relVel[2] * vCMag).to_value(u.km/u.s) )
            
            # Save the entire relativeEphem matrix.
            self.relPosArrayX = np.array(relPosArrayX) * ( 1 * u.km )
            self.relPosArrayY = np.array(relPosArrayY) * ( 1 * u.km )
            self.relPosArrayZ = np.array(relPosArrayZ) * (-1 * u.km )
            self.relVelArrayX = np.array(relVelArrayX) * ( 1 * u.km / u.s )
            self.relVelArrayY = np.array(relVelArrayY) * ( 1 * u.km / u.s )
            self.relVelArrayZ = np.array(relVelArrayZ) * (-1 * u.km / u.s )
            
            # To allow for chaining
            return self
        
        # Or if the orbit is a para/hyperbola...
        else:
            print('Error! Formations for non-closed orbits not supported!')
            return self
    
    def plot(self):
        
        """
        Function to plot the relative trajectory. You must run the propagate()
        method of the instance before the plot() method works.
        """
        
        # Check if the user has propagated the relative orbit
        if len(self.relPosArrayX) == 0 or len(self.relVelArrayX) == 0:
            print('The relative trajectories have not been propagated yet!')
            print('Have you forgotten to run the "propagate()" method?')
        
        # Else, proceed with the propagation!
        else:
            
            with quantity_support():
            
                # Initialise the matplotlib frame object.
                figMain = plt.figure()
                axOrbR = figMain.add_subplot(1,1,1, projection='3d')
                
                # Initialise the plot labels.
                axOrbR_label = 'Relative Position of Deputy'
                
                # Plot the relative orbit position in local VVLH frame.
                axOrbR.plot( self.relPosArrayZ, # X-track
                             self.relPosArrayY, # In-track
                             self.relPosArrayX, # Radial
                             label = axOrbR_label)
                
                # Set the relative orbit position axes labels.
                axOrbR.set_xlabel('Hill Frame Cross-Track Axis (km)')
                axOrbR.set_ylabel('Hill Frame In-Track Axis (km)')
                axOrbR.set_zlabel('Hill Frame Radial Axis (km)')
                
                # Get the current axes limits on relative orbit plots.
                axOrbR_axes_limits = [axOrbR.get_xlim()[0],
                                      axOrbR.get_xlim()[1],
                                      axOrbR.get_ylim()[0],
                                      axOrbR.get_ylim()[1],
                                      axOrbR.get_zlim()[0],
                                      axOrbR.get_zlim()[1]]
                
                # Using axOrbR_axes_limits above, we can find the minimum and
                # maximum axes limits and the span of values.
                axOrbR_axes_min  = min(axOrbR_axes_limits)
                axOrbR_axes_max  = max(axOrbR_axes_limits)
                axOrbR_axes_span = axOrbR_axes_max - axOrbR_axes_min
    
                # Scale all axes equally for relative orbit plots
                axOrbR.set_xlim( axOrbR_axes_min, axOrbR_axes_max )
                axOrbR.set_ylim( axOrbR_axes_min, axOrbR_axes_max )
                axOrbR.set_zlim( axOrbR_axes_min, axOrbR_axes_max )
                
                # It is important that the XYZ axes in the VVLH (relative 
                # orbit) frame is scaled the same, else it is difficult to 
                # interpret  the relative separations on different scales.
                
                # Plot the chief satellite as a tri-axial quiver in VVLH frame.
                axOrbR0 = axOrbR_axes_span * 0.1
                axOrbR.quiver( 0,0,0,1,0,0, length = axOrbR0,
                               color = 'r', arrow_length_ratio=0.3 )
                axOrbR.quiver( 0,0,0,0,1,0, length = axOrbR0,
                               color = 'r', arrow_length_ratio=0.3 )
                axOrbR.quiver( 0,0,0,0,0,1, length = axOrbR0,
                               color = 'r', arrow_length_ratio=0.3 )
                
                # Set tight layout.
                plt.tight_layout()
                plt.show()
    
    def plot_v(self):
        
        """
        Function to plot the relative velocity. You must run the propagate()
        method of the instance before the plot_v() method works.
        """
        
        # Check if the user has propagated the relative orbit
        if len(self.relPosArrayX) == 0 or len(self.relVelArrayX) == 0:
            print('The relative trajectories have not been propagated yet!')
            print('Have you forgotten to run the "propagate()" method?')
        
        # Else, proceed with the propagation!
        else:
            
            with quantity_support():
            
                # Initialise the matplotlib frame object.
                figMain = plt.figure()
                axOrbV = figMain.add_subplot(1,1,1, projection='3d')
                
                # Initialise the plot labels.
                axOrbV_label = 'Relative Velocity of Deputy'
                
                # Plot the relative orbit in local VVLH frame.
                axOrbV.plot( self.relVelArrayZ, # X-track
                             self.relVelArrayY, # In-track
                             self.relVelArrayX, # Radial
                             label = axOrbV_label)
                
                # Set the relative orbit position axes labels.
                axOrbV.set_xlabel('Hill Frame Cross-Track Velocity (km/s)')
                axOrbV.set_ylabel('Hill Frame In-Track Velocity (km/s)')
                axOrbV.set_zlabel('Hill Frame Radial Velocity (km/s)')
                
                # Get the current axes limits on relative orbit plots.
                axOrbV_axes_limits = [axOrbV.get_xlim()[0],
                                      axOrbV.get_xlim()[1],
                                      axOrbV.get_ylim()[0],
                                      axOrbV.get_ylim()[1],
                                      axOrbV.get_zlim()[0],
                                      axOrbV.get_zlim()[1]]
                
                # Using axOrbV_axes_limits above, we can find the minimum and
                # maximum axes limits and the span of values.
                axOrbV_axes_min  = min(axOrbV_axes_limits)
                axOrbV_axes_max  = max(axOrbV_axes_limits)
                axOrbV_axes_span = axOrbV_axes_max - axOrbV_axes_min
    
                # Scale all axes equally for relative orbit plots
                axOrbV.set_xlim( axOrbV_axes_min, axOrbV_axes_max )
                axOrbV.set_ylim( axOrbV_axes_min, axOrbV_axes_max )
                axOrbV.set_zlim( axOrbV_axes_min, axOrbV_axes_max )
                
                # Plot the chief satellite in-track velocity quiver in VVLH.
                axOrbV0 = axOrbV_axes_span * 0.1
                axOrbV.quiver( 0,0,0,1,0,0, length = axOrbV0,
                               color = 'r', arrow_length_ratio=0.3 )
                axOrbV.quiver( 0,0,0,0,1,0, length = axOrbV0,
                               color = 'r', arrow_length_ratio=0.3 )
                axOrbV.quiver( 0,0,0,0,0,1, length = axOrbV0,
                               color = 'r', arrow_length_ratio=0.3 )
                
                # Set tight layout.
                plt.tight_layout()
                plt.show()
            
##############################################################################
##############################################################################
##                                                                          ##
##  END OF THE RELATIVE ORBITS CLASS DEFINITION. CLASS INSTANTIATION BELOW. ##
##                                                                          ##
##############################################################################
##############################################################################

# Initialize an example Satellite A as the chief spacecraft.
satC = Orbit.from_classical( attractor = Earth,
                             a         = 6918.140 * u.km,
                             ecc       = 1e-6     * u.one,
                             inc       = 10.0     * u.deg,
                             raan      = 70.0     * u.deg,
                             argp      = 90.0     * u.deg,
                             nu        = 1.65     * u.deg)

# Initialize an example Satellite B as the deputy spacecraft.
satD = Orbit.from_classical( attractor = Earth,
                             a         = 6918.140 * u.km,
                             ecc       = 0.012    * u.one,
                             inc       = 11.4     * u.deg,
                             raan      = 72.35    * u.deg,
                             argp      = 135.0    * u.deg,
                             nu        = -46.5725 * u.deg)

# Instantiate the relative orbits object.
relativeSat = RelativeOrb( satC, satD )

# Propagate the relative orbit
relativeSat.propagate()

# Plot the relative trajectory in local Hill Frame (VVLH) of the chief.
relativeSat.plot()


