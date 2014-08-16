# -*- coding: utf-8 -*-
# _synthesizeNTF0.py
# Module providing the synthesizeNTF function
# Copyright 2013 Giuseppe Venturini
# This file is distributed with python-deltasigma.
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
# The following code has been (slightly) modified from pydsm, its original
# copyright notice follows:
#
# Copyright (c) 2012, Sergio Callegari
# All rights reserved.
#
# The py dsm code was ported from the MATLAB Delta Sigma toolbox
# Copyright (c) 2009, Richard Schreier
#
# The three software follow the same license, known as the 2-clause BSD.
# See the LICENSE file for details.

"""
Synthesize a noise transfer function (NTF) for a delta-sigma modulator without
optimizing the result.
"""

import numpy as np
from warnings import warn
from ._evalTF import evalTF
from ._utils import cplxpair
from ._ds_optzeros import ds_optzeros
from ._config import itn_limit

def synthesizeNTF0(order, osr, opt, H_inf, f0):
    """Synthesize a noise transfer function for a delta-sigma modulator.

    ::warn This function is not meant to be used directly, instead, set
    optimize_NTF to False and call synthesizeNTF(...).

    Parameters
    ----------
    order : int,
        the order of the modulator
    osr : float,
        the oversamping ratio
    opt : int or list of floats
        flag for optimized zeros

        * 0 -> not optimized,
        * 1 -> optimized,
        * 2 -> optimized with at least one zero at band-center,
        * [z] -> zero locations in complex form

    H_inf : real
        max allowed peak value of the NTF.
    f0 : real
        center frequency for BP modulators, or 0 for LP modulators.
        Note that 1 corresponds to the sampling frequency, so that 0.5 is the
        maximum value. Value 0 specifies an LP modulator.

    Returns
    -------
    ntf : tuple
        noise transfer function in zpk form.

    Raises
    ------
    ValueError
        'Cannot synthesize NTF zeros' if an empty list is supplied for the zeros.

        'Order must be even for a bandpass modulator.' if the order is
        incompatible with the modulator type.

        'The opt vector must be of length xxx' if opt is used to explicitly
        pass the NTF zeros and these are in the wrong number.

    Warns
    -----
        'Unable to achieve specified H_inf ...' if the desired H_inf
        cannot be achieved.

        'Iteration limit exceeded' if the routine converges too slowly.

    Notes
    -----
    This is actually a wrapper function which calls the appropriate version
    of synthesizeNTF, based on the module control flag `optimize_NTF` which
    determines whether to use optimization tools.

    Parameter H_inf is used to enforce the Lee stability criterion.

    See also:
       clans()   "Closed-loop analysis of noise-shaper." An alternative
                 method for selecting NTFs based on the 1-norm of the
                 impulse response of the NTF

       synthesizeChebyshevNTF()    Select a type-2 highpass Chebyshev NTF.
                 This function does a better job than synthesizeNTF if osr
                 or H_inf is low.
    """

    # Determine the zeros.
    if f0 != 0:
        # Bandpass design-- halve the order temporarily.
        order = order/2
        dw = np.pi/(2*osr)
    else:
        dw = np.pi/osr

    if np.isscalar(opt):
        # opt is a number
        if opt == 0:
            z = np.zeros(order)
        else:
            z = dw*ds_optzeros(order, opt)
        if z.size == 0:
            raise ValueError('Cannot synthesize NTF zeros')
        if f0 != 0:
            # Bandpass design-- shift and replicate the zeros.
            order = order*2
            z = z + 2*np.pi*f0
            z = np.vstack((z, -z)).transpose().flatten()
        z = np.exp(1j*z)
    else:
        z = opt

    p = np.zeros(order)
    k = 1
    fprev = 0

    if f0 == 0:
        # Lowpass design
        HinfLimit = 2**order
        # !!! The limit is actually lower for opt=1 and low OSR
        if H_inf >= HinfLimit:
            warn('Unable to achieve specified H_inf.\n'
                 'Setting all NTF poles to zero.')
            p = np.zeros(order)
        else:
            x = 0.3**(order-1)   # starting guess
            for itn in range(1, itn_limit + 1):
                me2 = -0.5*(x**(2./order))
                w = (2*np.arange(1, order + 1) + 1)*np.pi/order
                mb2 = 1 + me2*np.exp(1j*w)
                p = mb2 - np.sqrt(mb2**2 - 1)
                # Reflect poles to be inside the unit circle
                out = abs(p) > 1
                p[out] = 1./p[out]
                # The following is not exactly what delsig does.
                # We do not have an identical cplxpair
                p = cplxpair(p)
                f = np.real(evalTF((z, p, k), -1))-H_inf
                if itn == 1:
                    delta_x = -f/100.
                else:
                    delta_x = -f*delta_x/(f - fprev)
                xplus = x + delta_x
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
                if itn == itn_limit:
                    warn('Iteration limit exceeded.')
    else:
        # Bandpass design
        x = 0.3**(order/2-1)   # starting guess (not very good for f0~0)
        if f0 > 0.25:
            z_inf = 1.
        else:
            z_inf = -1.
        c2pif0 = np.cos(2*np.pi*f0)
        for itn in range(1, itn_limit+1):
            e2 = 0.5*x**(2./order)
            w = (2*np.arange(order)+1)*np.pi/order
            mb2 = c2pif0 + e2*np.exp(1j*w)
            p = mb2 - np.sqrt(mb2**2-1)
            # Reflect poles to be inside the unit circle
            out = abs(p)>1
            p[out] = 1/p[out]
            # The following is not exactly what delsig does.
            p = cplxpair(p)
            f = np.real(evalTF((z, p, k), z_inf)) - H_inf
            if itn == 1:
                delta_x = -f/100
            else:
                delta_x = -f*delta_x/(f - fprev)
            xplus = x + delta_x
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
            if itn == itn_limit:
                warn('Iteration limit exceeded.')

    z = cplxpair(z)
    return z, p, k

