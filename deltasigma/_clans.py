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

"""Module providing the clans() optimization function.
"""

from __future__ import division, print_function
import numpy as np
from scipy.optimize import minimize
from scipy.signal import dimpulse

from ._synthesizeNTF import synthesizeNTF
from ._utils import cplxpair
from ._evalTF import evalTF
from ._dsclansNTF import dsclansNTF

# We ported the version that was originally designed for the MATLAB
# Optimization Toolbox version >= 6

def clans(order=4, OSR=64, Q=5, rmax=0.95, opt=0):
    """Optimal NTF design for a multi-bit modulator.
    
    Synthesize a noise transfer function (NTF) for a lowpass delta-sigma 
    modulator using the CLANS methodology.

    CLANS stands for Closed-Loop Analysis of Noise-Shapers,
    and it was originally developed by J.G. Kenney and L.R. Carley [1]_.

    .. [1] J. G. Kenney and L. R. Carley, “Design of multibit
        noise-shaping data converters,” Analog Integrated Circuits Signal
        Processing Journal, vol. 3, pp. 259-272, 1993.

    **Parameters:**

    order : int
        The order of the NTF.

    OSR : int
        The oversampling ratio.

    Q : int
        The maximum number of quantization levels used by the fed-back quantization
        noise. (Mathematically, :math:`Q = \\|h\\|_1 - 1`, i.e. the sum of the 
        absolute values of the impulse response samples minus one is the 
        maximum instantaneous noise gain.)

    rmax : float
        The maximum radius for the NTF poles.

    opt : int
        A flag used to request optimized NTF zeros. 
        
    * `opt=0` puts all NTF zeros at band center (DC for lowpass modulators).
    * `opt=1` optimizes the NTF zeros.
    * For even-order modulators, `opt=2` puts two zeros at band-center, 
      but optimizes the rest.

    **Returns**

    ntf : tuple
        The modulator NTF, given in ZPK (zero-pole-gain) form.

    **Example:**

    Fifth-order lowpass modulator; (time-domain) noise gain of 5, 
    zeros optimized for OSR = 32.::

        H = clans(5, 32, 5, .95, 1)
        pretty_lti(H)

    Returns::

              (z -1) (z^2 -1.997z +1) (z^2 -1.992z +0.9999)      
        ---------------------------------------------------------
         (z -0.4184) (z^2 -1.305z +0.5713) (z^2 -0.978z +0.2686) 

    ``H`` can be plotted through :func:`DocumentNTF`:

    .. plot::

        from deltasigma import DocumentNTF, clans
        # Fifth-order lowpass modulator; (time-domain) noise gain of 5, 
        # zeros optimized for OSR = 32.
        H = clans(5, 32, 5, .95, 1)
        DocumentNTF(H)

    """
    # Create the initial guess
    Hz, poles, _ = synthesizeNTF(order, OSR, opt, 1 + Q, 0)
    x = np.zeros((order, ))
    odd = order % 2
    poles = cplxpair(poles)
    poles = poles[::-1]
    if odd == 1:
        z = poles[0]/rmax
        if (np.abs(z) > 1).any(): #project poles outside rmax onto the circle
            z = z/np.abs(z)
        s = (z - 1)/(z + 1)
        x[0] = np.real_if_close(np.sqrt(-s))
    for i in range(odd, order, 2):
        z = poles[i:i + 2]/rmax
        if np.any(np.abs(z) > 1): #project poles outside rmax onto the circle
            z = z/np.abs(z)
        s = (z - 1)/(z + 1)
        coeffs = np.poly(s)
        wn = np.sqrt(coeffs[2])
        zeta = coeffs[1]/(2 * wn)
        x[i] = np.real_if_close(np.sqrt(zeta))
        x[i + 1] = np.real_if_close(np.sqrt(wn))
    #options=optimset('TolX',1e-06,'TolFun',1e-06,'TolCon',1e-06,'MaxIter',1000)
    #fobj = lambda x: dsclansObj6a(x, order, OSR, Q, rmax, Hz)
    #fconstr = lambda x: dsclansObj6b(x, order, OSR, Q, rmax, Hz)
    res = minimize(dsclansObja, x, args=(order, OSR, Q, rmax, Hz), 
             method='slsqp', constraints={'type':'ineq', 'fun':dsclansObjb, 
             'args':(order, OSR, Q, rmax, Hz)})
    NTF = dsclansNTF(res['x'], order, rmax, Hz)
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
    H = (H[0], H[1], H[2], 1.)
    # dimpulse(H, n=100)[y is 0][output 0]
    g = np.sum(np.abs(dimpulse(H, t=np.arange(100))[1][0])) - 1 - Q
    return -g

