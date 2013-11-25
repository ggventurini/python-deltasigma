# -*- coding: utf-8 -*-
# _clans.py
# Closed-Loop Analysis of Noise-Shapers module
# Copyright 2013 Giuseppe Venturini
# This file is part of python-deltasigma.
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's 
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.

from __future__ import division
from warnings import warn
import numpy as np
from scipy.optimize import minimize
from scipy.signal import impulse

from ._synthesizeNTF import synthesizeNTF
from ._utils import cplxpair
from ._evalTF import evalTF
from ._dsclansNTF import dsclansNTF

"""Module providing the clans() optimization function."""

# We ported the version that was originally designed for the MATLAB
# Optimization Toolbox version >= 6

def clans(order=4, OSR=64, Q=5, rmax=0.95, opt=0):
    """Optimal NTF design for a multi-bit modulator.
    
    CLANS stands for Closed-Loop Analysis of Noise-Shapers,
    and it was originally developed by J.G. Kenney and L.R. Carley.
    """
    warn("Warning untested function.")
    # Create the initial guess
    Hz, poles, _ = synthesizeNTF(order, OSR, opt, 1 + Q, 0)
    x = np.zeros((1, order))
    odd = order % 2
    poles = cplxpair(poles)
    poles = poles[::-1]
    if odd == 1:
        z = poles[0]/rmax
        if (np.abs(z) > 1).any(): #project poles outside rmax onto the circle
            z = z/abs(z)
        s = (z - 1)/(z + 1)
        x[0, 0] = np.sqrt(-s)
    for i in range(odd, order, 2):
        z = poles[i:i + 1] / rmax
        if (np.abs(z) > 1).any(): #project poles outside rmax onto the circle
            z = z/np.abs(z)
        s = (z - 1)/(z + 1)
        coeffs = np.poly(s)
        wn = np.sqrt(coeffs[2])
        zeta = coeffs[1]/(2 * wn)
        x[i] = np.sqrt(zeta)
        x[i + 1] = np.sqrt(wn)
    #options=optimset('TolX',1e-06,'TolFun',1e-06,'TolCon',1e-06,'MaxIter',1000)
    #fobj = lambda x: dsclansObj6a(x, order, OSR, Q, rmax, Hz)
    #fconstr = lambda x: dsclansObj6b(x, order, OSR, Q, rmax, Hz)
    minimize(dsclansObja, x, args=(order, OSR, Q, rmax, Hz), 
             cons={'type':'ineq', 'fun':dsclansObjb, 
             'args':(order, OSR, Q, rmax, Hz)})
    NTF = dsclansNTF(x, order, rmax, Hz)
    return NTF

def dsclansObja(x, order, OSR, Q, rmax, Hz):
    """Objective function for clans; Optimization Toolbox version >= 6
    f is the magnitude of H at the band-edge
    """
    H = dsclansNTF(x, order, rmax, Hz)
    f = np.abs(evalTF(H, np.exp(1j*np.pi/OSR)))
    return f

def dsclansObjb(x, order, OSR, Q, rmax, Hz):
    """Constraint function for clans; 
    g =||h||_1 - Q
    """
    H = dsclansNTF(x, order, rmax, Hz)
    g = np.sum(np.abs(impulse(H, 100))) - 1 - Q
    return g
