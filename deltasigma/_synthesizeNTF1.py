# -*- coding: utf-8 -*-
# _synthesizeNTF1.py
# Module providing the synthesizeNTF1 function
# This file is distributed with python-deltasigma.
# Copyright 2013 Giuseppe Venturini
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.
#
# The following code has been taken with little modifications from pydsm,
# its original copyright notice follows:
#
# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# The code was ported from the MATLAB Delta Sigma toolbox, which is
# Copyright (c) 2009, Richard Schreier
#
# The three software follow the same license, known as the 2-clause BSD.
# See the LICENSE file for details.

"""
Module providing the synthesizeNTF1() function.
"""

# The following code is
# Copyright (c) 2012, Sergio Callegari
# All rights reserved.

# Portions of code ported from the DELSIG toolbox
# Copyright (c) 2009, Richard Schreier

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from warnings import warn
from ._evalTF import evalTF
from ._utils import cplxpair
from ._ds_optzeros import ds_optzeros
from ._ds_synNTFobj1 import ds_synNTFobj1
from ._padl import padl
from ._constants import eps

def synthesizeNTF1(order, osr, opt, H_inf, f0):
    """
    Synthesize a noise transfer function (NTF) for a delta-sigma modulator optimizing the result.
    """
    # Determine the zeros.
    if f0 != 0:
        # Bandpass design-- halve the order temporarily.
        order = order/2
        dw = np.pi/(2*osr)
    else:
        dw = np.pi/osr

    if opt.ndim == 0:
        # opt is a number
        if opt == 0:
            z = np.zeros(order)
        else:
            z = dw*ds_optzeros(order, 1 + np.fmod(opt-1, 2))
        if z.size == 0:
            raise ValueError('Cannot synthesize NTF zeros')
        if f0 != 0:
            # Bandpass design-- shift and replicate the zeros.
            order = order*2
            z = np.sort(z) + 2*np.pi*f0
            z = np.vstack((z,-z)).transpose().flatten()
        z = np.exp(1j*z)
    else:
        z = opt

    zp = z[np.angle(z) > 0]
    x0 = (np.angle(zp)-2*np.pi*f0) * osr / np.pi
    if opt.size == 1 and opt == 4 and f0 != 0:
        # Do not optimize the zeros at f0
        x0 = np.delete(x0, np.nonzero(abs(x0) < eps))

    p = np.zeros(order)
    k = 1
    Hinf_itn_limit = 100
    fprev = 0

    opt_iteration = 5   # Max number of zero-optimizing/Hinf iterations
    while opt_iteration > 0:
        # Iteratively determine the poles by finding the value of the x
        # parameter which results in the desired H_inf
        ftol = 1e-10
        if f0 > 0.25:
            z_inf = 1
        else:
            z_inf = -1
        if f0 == 0:
            # Lowpass design
            HinfLimit = 2**order
            # !!! The limit is actually lower for opt=1 and low OSR
            if H_inf >= HinfLimit:
                warn('Unable to achieve specified Hinf.\n'
                    'Setting all NTF poles to zero.')
                p = np.zeros(order)
            else:
                x = 0.3**(order-1)   # starting guess
                for itn in range(1, Hinf_itn_limit+1):
                    me2 = -0.5*(x**(2./order))
                    w = (2*np.arange(1,order+1)+1)*np.pi/order
                    mb2 = 1+me2*np.exp(1j*w)
                    p = mb2 - np.sqrt(mb2**2-1)
                    # Reflect poles to be inside the unit circle
                    out = abs(p)>1
                    p[out] = 1/p[out]
                    # The following is not exactly what delsig does.
                    # We do not have an identical cplxpair
                    p = cplxpair(p)
                    f = np.real(evalTF((z, p, k), z_inf))-H_inf
                    if itn == 1:
                        delta_x = -f/100
                    else:
                        delta_x = -f*delta_x/(f-fprev)

                    xplus = x+delta_x
                    if xplus > 0:
                        x = xplus
                    else:
                        x = x*0.1
                    fprev = f
                    if abs(f) < 1e-10 or abs(delta_x) < 1e-10:
                        break
                    if x > 1e6:
                        warn('Unable to achieve specified Hinf.\n'
                             'Setting all NTF poles to zero.')
                        p = np.zeros(order)
                        break
                    if itn == Hinf_itn_limit:
                        warn('Danger! Iteration limit exceeded.')
        else:
            # Bandpass design
            x = 0.3**(order/2-1)   # starting guess (not very good for f0~0)
            if f0 > 0.25:
                z_inf = 1.
            else:
                z_inf = -1.
            c2pif0 = np.cos(2*np.pi*f0)
            for itn in range(1, Hinf_itn_limit+1):
                e2 = 0.5*x**(2./order)
                w = (2*np.arange(order)+1)*np.pi/order
                mb2 = c2pif0 + e2*np.exp(1j*w)
                p = mb2 - np.sqrt(mb2**2-1)
                # Reflect poles to be inside the unit circle
                out = abs(p)>1
                p[out] = 1/p[out]
                # The following is not exactly what delsig does.
                p = cplxpair(p)
                f = np.real(evalTF((z, p, k), z_inf))-H_inf
                if itn == 1:
                    delta_x = -f/100
                else:
                    delta_x = -f*delta_x/(f-fprev)
                xplus = x+delta_x
                if xplus > 0:
                    x = xplus
                else:
                    x = x*0.1
                fprev = f
                if abs(f) < 1e-10 or abs(delta_x) < 1e-10:
                    break
                if x > 1e6:
                    warn('Unable to achieve specified Hinf.\n'
                        'Setting all NTF poles to zero.')
                    p = np.zeros(order)
                    break
                if itn == Hinf_itn_limit:
                    warn('Danger! Iteration limit exceeded.')

        # ---- Zero optimization part
        if (opt.size == 1 and opt < 3) or opt.size > 1:
            # Do not optimize the zeros
            opt_iteration = 0
        else:
            if f0 == 0:
                ub = np.ones(x0.size)
                lb = np.zeros(x0.size)
            else:
                ub = 0.5*np.ones(x0.size)
                lb = -ub
            # options = optimset('TolX',0.001, 'TolFun',0.01, 'MaxIter',100 );
            # options = optimset(options,'LargeScale','off');
            # options = optimset(options,'Display','off');
            # %options = optimset(options,'Display','iter');
            opt_result = fmin_l_bfgs_b(ds_synNTFobj1, x0, args=(p, osr, f0),
                                      approx_grad=True, bounds=list(zip(lb,ub)))
            x=opt_result[0]
            x0 = x
            z = np.exp(2j*np.pi*(f0+0.5/osr*x))
            if f0 > 0:
                z = padl(z, len(p)/2, np.exp(2j*np.pi*f0))
            z = np.concatenate((z, z.conj()), axis=1)
            if f0 == 0:
                z = padl(z, len(p), 1)
            if  np.abs(np.real(evalTF((z, p, k), z_inf)) - H_inf ) < ftol:
                opt_iteration = 0
            else:
                opt_iteration = opt_iteration - 1
    z = cplxpair(z)
    return (z, p, k)
