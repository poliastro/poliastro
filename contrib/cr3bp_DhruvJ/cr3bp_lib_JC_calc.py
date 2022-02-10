# -*- coding: utf-8 -*-
"""
Created on 7 Jan 22 12:24:03 2021

@author: Dhruv Jain, Multi-Body Dynamics Research Group, MSAAE Purdue University
        dhruvj9922@gmail.com

Objectve: This function has the two following features:  
    1. Calculate the position [nd] of 5 libration points for a system in CR3BP
    2. Compute Jacobi Constant [nd] 
"""
import numpy as np

def lib_pt_loc(mu):
    """ Computes Non-Dimensionalized Libration Points Location for P1-P2 system
    Args: 
    mu: float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    Returns: 
    lib_loc: numpy ndarray (5x3)   
        5 Libration Points, [nd]
    """

    lib_loc = np.zeros((5,3))
    lib_loc[3,:] = [0.5-mu, 3**0.5/2, 0] # L4
    lib_loc[4,:] = [0.5-mu, -3**0.5/2, 0] # L5

    # 5th degree polynomial of L1, L2 and L3 
    f_lib = np.array([[1, mu-3, 3-2*mu, -mu, 2*mu, -mu], 
                      [1, 3-mu, 3-2*mu, -mu, -2*mu, -mu],
                      [1, 2+mu, 1+2*mu, mu-1, 2*mu-2, -1+mu]])
    
    fd_lib = np.array([[0, 5, 4*(mu-3), 3*(3-2*mu), 2*-mu, 2*mu],
                      [0, 5, 4*(3-mu), 3*(3-2*mu), 2*-mu, -2*mu],
                      [0, 5, 4*(2+mu), 3*(1+2*mu), 2*(mu-1), 2*mu-2]]) 
    
    ic = np.array([0.9, 1.1, -1])
    
    tolerance = 1e-10
    
    for i in range(3):
        val = np.vander([ic[i]],6)
        h = np.dot(val,f_lib[i,:])/np.dot(val,fd_lib[i,:])
        while abs(h) >=tolerance:
            val = np.vander([ic[i]],6)
            h = np.dot(val,f_lib[i,:])/np.dot(val,fd_lib[i,:])          
            lib_loc[i,0] = ic[i]-h
            
            ic[i] = lib_loc[i,0]            
        
        if i == 0:
            lib_loc[i,0] = 1-mu-lib_loc[i,0]  
        elif i == 1:
            lib_loc[i,0] = 1-mu+lib_loc[i,0]  
        elif i == 2:
            lib_loc[i,0] = -mu-lib_loc[i,0]  
            
            
    return lib_loc

def JC(mu,r,v):
    """ Computes Jacobi Constant/Jacobi Integral [nd], CR3BP integral constant
    It can handle complex inputs
    Args:
        mu: float, M2/(M1+M2)
            M1 and M2 are mass of Primary Bodies and M2<M1
        r: numpy ndarray, 3x1
           position vector (x,y,z), [nd]
        v: numpy ndarray, 3x1
           veloctiy vector (vx, vy, vz), [nd]
    Returns:
        JC: Jacobi Constant, [nd]
        
    """
    d1 = ((r[0]+mu)**2 + r[1]**2 + r[2]**2)**0.5 # d1 = d
    d2 = ((r[0]-1+mu)**2 + r[1]**2 + r[2]**2)**0.5 # d2 = r 
    JC = r[0]**2 + r[1]**2 + 2*(1-mu)/d1 + 2*mu/d2 - np.sqrt(v[0]**2+v[1]**2+v[2]**2)**2 # linalg.norm expanded for complex v
    
    return JC
    
    
    
    
    
    
    
    
    
    
    
    
    
    
