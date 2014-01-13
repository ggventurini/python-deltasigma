# -*- coding: utf-8 -*-
# _synthesizeChebyshevNTF.py
# Module providing the synthesizeChebyshevNTF function
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

"""Module providing the synthesizeChebyshevNTF() function
"""

from __future__ import division
from warnings import warn

import numpy as np

from scipy.signal import cheby2

from ._ds_f1f2 import ds_f1f2

def synthesizeChebyshevNTF(order=3, OSR=64, opt=0, H_inf=1.5, f0=0.):
    """Synthesize a noise transfer function for a delta-sigma modulator.
    
    The NTF is a type-2 highpass Chebyshev function.

    **Parameters:**

    order : int, optional
        order of the modulator, defaults to 3

    OSR : int, optional
        oversampling ratio, defaults to 64

    opt : int, optional
        ignored value, for consistency with ::func:synthesizeNTF

    H_inf : float, optional
        maximum NTF gain, defaults to 1.5

    f0 : float, optional
        center frequency (1->fs), defaults to 0.

    **Returns:**

    z, p, k : tuple 
        a zpk tuple containing the zeros and poles of the NTF.

    **Warns:**

    * If a non-zero value is passed for ``opt``.

    **Raises:**

    * ValueError: Order must be even for a bandpass modulator

    """
    if opt:
        warn("Got a non-zero 'opt' value. Not such optimization is " + \
             "available, opt is only meant to ease switching between " + \
             "synthesizeNTF and synthesizeChebyshevNTF.")
    if f0 != 0:
        if order % 2 != 0:
            raise ValueError('Order must be even for a bandpass modulator.')
        else:
            f1, f2 = ds_f1f2(OSR, f0)
            f1f2 = np.array([f1, f2])
    x_min = 0
    x_max = 300
    dx_max = 10
    ftol = 1e-06
    xtol = 1e-06
    x = 60
    itn_limit = 10
    converged = False
    for itn in range(itn_limit):
        if f0 == 0:
            z, p, k = cheby2(order, x, 1./OSR, btype='high', output='zpk')
        else:
            z, p, k = cheby2(order/2., x, 2.*f1f2, btype='stop', output='zpk')
        f = 1./k - H_inf
        if f > 0:
            x_max = x
        else:
            x_min = x
        if itn == 0:
            dx = -dx_max*np.sign(f)
        else:
            df = f - f_p
            if abs(df) < ftol:
                converged = True
                break
            dx = -f*dx/df
        if converged:
            break
        x_p = x
        f_p = f
        x = max(x_min, min(x + dx, x_max))
        dx = x - x_p
        if abs(dx) < xtol:
            break
    ntf = (z, p, 1)
    return ntf

def test_synthesizeChebyshevNTF():
    """Test function for synthesizeChebyshevNTF()"""
    from ._utils import cplxpair
    z, p, k = synthesizeChebyshevNTF()
    zref = [1., .9991 + 0.0425j, .9991 - 0.0425j]
    pref = [.6609, .7686 + .2858j, .7686 - .2858j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    z, p, k = synthesizeChebyshevNTF(order=4, OSR=32, opt=0, H_inf=1.5, f0=.33)
    zref = [-.4513 + .8924j, -.4513 - .8924j, -.5122 + 0.8589j, -.5122 - 0.8589j]
    pref = [-.2249 + .7665j, -.2249 - .7665j, -.5506 + .6314j, -.5506 - .6314j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
