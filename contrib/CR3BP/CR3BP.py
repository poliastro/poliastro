import numpy as np
from numba import njit as jit

from poliastro._math.ivp import DOP853, solve_ivp


@jit
def getChar_CR3BP(k1,k2,r12):
    """Characteristic values for Circular Restricted Three Body Problem (CR3BP).
    All parameters are in non-dimensional form. 
    
    Parameters
    ----------
    k1 : float
        Primary body gravitational parameter.
    k2 : float
        Secondary body gravitational parameter.
    r12 : float
        mean distance between two bodies
    Returns
    ----------
    mu : float
        CR3BP mass ratio.
    kstr : float
        Characterisitc gravitational parameter.
    lstr : float
        Characterisitc length.
    tstr : float
        Characterisitc time.
    vstr : float
        Characterisitc velocity.
    nstr : float
        Characterisitc angular velocity.
    """
    
    # characteristic gravitational parameter
    kstr = k1 + k2
    
    # characteristic mass
    # mstr = kstr/G
    
    # characteristic length
    lstr = r12
    
    # characteristic time
    tstr = (lstr**3/kstr) ** 0.5
    
    # characteristic velocity
    vstr = lstr/tstr
    
    # characteristic angular velocity
    nstr = (kstr/r12**3) ** 0.5
    
    # CR3BP mass ratio
    mu = k2/kstr
    
    return mu, kstr, lstr, tstr, vstr, nstr


@jit
def getJacobian_CR3BP(u_, mu):
    """Compute state matrix A of Circular Restricted Three Body Problem (CR3BP).
    All parameters are in non-dimensional form.
    
    Parameters
    ----------
    u_ : numpy.ndarray
        Six component state vector [rx, ry, rz, vx, vy, vz] (nd).
    mu : float
        CR3BP mass ratio (m2/(m1+m2))
    """
    rx, ry, rz, vx, vy, vz = u_
    
    # distance to primary body
    r13 = ( (rx+mu)**2 + ry**2 + rz**2 )**0.5;
    # distance to secondary body
    r23 = ( (rx-(1-mu))**2 + ry**2 + rz**2 )**0.5;
    
    # computing velocity square
    v2 = vx**2 + vy**2 + vz**2
    
    # compute Jacobi constant for trajectory
    C =  (rx**2+ry**2) + 2*(1-mu)/r13 + 2*mu/r23 - v2;
    
    return C

@jit
def getUdiff_CR3BP(r_, mu):
    """Compute Hetian of pseudo-potential of Circular Restricted Three Body 
    Problem (CR3BP).
    All parameters are in non-dimensional form.
    
    Parameters
    ----------
    r_ : numpy.ndarray
        Three component position vector [rx, ry, rz] (nd).
    mu : float
        CR3BP mass ratio (m2/(m1+m2))
    """
    
    # extracting components of position
    rx, ry, rz = r_
    
    # distance to primary body
    r13 = ( (rx+mu)**2 + ry**2 + rz**2 )**0.5;
    # distance to secondary body
    r23 = ( (rx-(1-mu))**2 + ry**2 + rz**2 )**0.5;
    
    # computing the double derivates with position
    Uxx = 1 - (1-mu)/r13**3 - mu/r23**3 
    + 3*(1-mu)*(rx+mu)**2/r13**5 + 3*mu*(rx+mu-1)**2/r23**5
    Uyy = 1 - (1-mu)/r13**3 - mu/r23**3 + 3*(1-mu)*ry**2/r13**5
    + 3*mu*ry**2/r23**5
    Uzz = -(1-mu)/r13**3 - mu/r23**3 + 3*(1-mu)*rz**2/r13**5 
    + 3*mu*rz**2/r23**5
    Uxy = 3*ry*(1-mu)*(rx+mu)/r13**5 + 3*ry*mu*(rx-(1-mu))/r23**5;
    Uxz = 3*rz*(1-mu)*(rx+mu)/r13**5 + 3*rz*mu*(rx-(1-mu))/r23**5;
    Uyz = 3*ry*rz*(1-mu)/r13**5 + 3*ry*rz*mu/r23**5;
    
    # exploiting the symmetry
    Uyx = Uxy
    Uzx = Uxz
    Uzy = Uyz
    
    # final Hetian matrix
    Udiff = np.array( [ [Uxx,Uxy,Uxz], [Uyx,Uyy,Uyz], [Uzx,Uzy,Uzz] ] )
    
    return Udiff


@jit
def getA_CR3BP(r_, mu):
    """Compute state matrix A of Circular Restricted Three Body Problem (CR3BP).
    All parameters are in non-dimensional form.
    
    Parameters
    ----------
    r_ : numpy.ndarray
        Three component position vector [rx, ry, rz] (nd).
    mu : float
        CR3BP mass ratio (m2/(m1+m2))
    """
    
    # Compute Hetian of psuedo-potential function(U)
    Udiff = getUdiff_CR3BP(r_,mu)
    
    # intialize A matrix 
    A = np.zeros((6,6))
    
    # adding upper left zeros
    A[0:3,0:3] = np.zeros((3,3))
    # adding upper right identity matrix
    A[0:3,3:6] = np.ones((3,3))
    # adding lower left hetian of pseudo-potential
    A[3:6,0:3] = Udiff
    # adding lower right Omega matrix
    A[3:6,3:6] = np.array([[0,2,0], [-2,0,0], [0,0,0]])
    
    return A


@jit
def func_CR3BP(t, u_, mu):
    """Differential equation for the initial value Circular Restricted Three
    Body Problem (CR3BP).
    All parameters are in non-dimensional form. 
    
    Parameters
    ----------
    t  : float
        Time (nd).
    u_ : numpy.ndarray
        Six component state vector [rx, ry, rz, vx, vy, vz] (nd).
    mu : float
        CR3BP mass ratio (m2/(m1+m2))
    """
    
    # extracting states
    rx, ry, rz, vx, vy, vz = u_
    
    # distance to primary body
    r13 = ( (rx+mu)**2 + ry**2 + rz**2 )**0.5;
    # distance to secondary body
    r23 = ( (rx-(1-mu))**2 + ry**2 + rz**2 )**0.5;
    
    # computing three-body dyamics
    r_dot = np.array([vx, vy, vz]);
    v_dot = np.array([
            rx + 2*vy - (1-mu)*(rx+mu)/(r13**3) - mu*(rx-1+mu)/(r23**3),
            ry - 2*vx - (1-mu)*ry/(r13**3) - mu*ry/(r23**3),
            -(1-mu)*rz/(r13**3) - mu*rz/(r23**3)
            ]);
    
    
    # state derivatives
    du = np.append(r_dot,v_dot)
    
    return du

@jit
def func_STM(t, u_, mu):
    
    # extract STM from dynamics vector 6:42 elements
    STM = u_[6:].reshape(6,6) # reshaped to a 6x6 matrix
    
    # get A matrix
    A = getA_CR3BP(u_[:3], mu)
    
    # compute STM derivatie
    STMdot = A @ STM
    
    # convert derivative matrix to a vector
    du = STMdot.reshape(1,36)
    
    return du


@jit
def func_CR3BP_STM(t, u_, mu):
    """Differential equation for the initial value Circular Restricted Three
    Body Problem (CR3BP) with State Transition Matrix Propagation.
    All parameters are in non-dimensional form. 
    
    Parameters
    ----------
    t  : float
        Time (nd).
    u_ : numpy.ndarray
        42 component vector [state + STM] (nd).
    mu : float
        CR3BP mass ratio (m2/(m1+m2))
    """
    
    # compute CR3BP state dynamics
    du_state = func_CR3BP(t,u_[0:6],mu)
    
    # compute CR3BP STM dynamics
    du_STM = func_STM(t,u_,mu)
    
    # full derivative vector
    du = np.append(du_state,du_STM)
    
    return du


def propagate(mu, r0, v0, tofs, rtol=1e-11, f=func_CR3BP):
    """Propagate an CR3BP orbit some time and return the result.
    
    Parameters
    ----------
    mu : float
        CR3BP mass ratio (nd)
    r0 : numpy.ndarray
        Position vector (nd).
    v0 : numpy.ndarray
        Velocity vector (nd).
    tofs : numpy.ndarray
        Array of times to propagate (nd).
    rtol : float, optional
        Maximum relative error permitted, defaults to 1e-11.
    f : function(t0, u, k), optional
        Objective function, default to natural CR3BP modell.
        
    Returns
    -------
    rr : numpy.ndarray
        Propagated position vectors (nd).
    vv : numpy.ndarray
        Propagated velocity vectors (nd).
    """

    u0 = np.append(r0,v0)

    result = solve_ivp(
        f,
        (0, max(tofs)),
        u0,
        args=(mu,),
        rtol=rtol,
        atol=1e-12,
        method=DOP853,
        dense_output=True
    )
    
    if not result.success:
        raise RuntimeError("Integration failed")


    rrs = []
    vvs = []
    for i in range(len(tofs)):
        t = tofs[i]
        y = result.sol(t)
        rrs.append(y[:3])
        vvs.append(y[3:])

    return rrs, vvs


def propagateSTM(mu, r0, v0, STM0, tofs, rtol=1e-11, f=func_CR3BP_STM):
    """Propagate an CR3BP orbit with STM some time and return the result.
    
    Parameters
    ----------
    mu : float
        CR3BP mass ratio (nd)
    r0 : numpy.ndarray
        Position vector (nd).
    v0 : numpy.ndarray
        Velocity vector (nd).
    STM0 : numpy.ndarray
        6x6 intial STM matrix (nd)
    tofs : numpy.ndarray
        Array of times to propagate (nd).
    rtol : float, optional
        Maximum relative error permitted, defaults to 1e-11.
    f : function(t0, u, k), optional
        Objective function, default to natural CR3BP modell.
        
    Returns
    -------
    rr : numpy.ndarray
        Propagated position vectors (nd).
    vv : numpy.ndarray
        Propagated velocity vectors (nd).
    STM : numpy.ndarray
        Propagated STM matrix (nd)
    """
    
    u0 = np.append(r0,v0)
    u0 = np.append(u0,STM0.reshape(1,36))

    result = solve_ivp(
        f,
        (0, max(tofs)),
        u0,
        args=(mu,),
        rtol=rtol,
        atol=1e-12,
        method=DOP853,
        dense_output=True
    )
    
    if not result.success:
        raise RuntimeError("Integration failed")


    rrs = []
    vvs = []
    STMs = []
    for i in range(len(tofs)):
        t = tofs[i]
        y = result.sol(t)
        rrs.append(y[:3])
        vvs.append(y[3:6])
        STMs.append(y[6:].reshape(6,6))

    return rrs, vvs, STMs