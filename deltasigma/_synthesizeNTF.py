# -*- coding: utf-8 -*-
# _synthesizeNTF.py
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
# The three software follow the same license, known as the 3-clause BSD. 
# See the LICENSE file for details.

"""
Module providing the main NTF synthesis function.
"""

import numpy as np
from warnings import warn
from ._config import optimize_NTF
from ._synthesizeNTF0 import synthesizeNTF0
from ._synthesizeNTF1 import synthesizeNTF1

def synthesizeNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0.0):
    """Synthesize a noise transfer function for a delta-sigma modulator.

    **Parameters:**

    order : int, optional
        the order of the modulator, defaults to 3
    osr : float, optional
        the oversamping ratio, defaults to 64
    opt : int or list of floats, optional
        flag for optimized zeros, defaults to 0

        * 0 -> not optimized,
        * 1 -> optimized,
        * 2 -> optimized with at least one zero at band-center,
        * 3 -> optimized zeros (with optimizer)
        * 4 -> same as 3, but with at least one zero at band-center
        * [z] -> zero locations in complex form

    H_inf : real, optional
        max allowed peak value of the NTF. Defaults to 1.5
    f0 : real, optional
        center frequency for BP modulators, or 0 for LP modulators.
        Defaults to 0.
        1 corresponds to the sampling frequency, so that 0.5 is the
        maximum value. Value 0 specifies an LP modulator.

    **Returns:**

    ntf : tuple
        noise transfer function in zpk form.

    **Raises:**

    ValueError
        * 'Error. f0 must be less than 0.5' if f0 is out of range

        * 'Order must be even for a bandpass modulator.' if the order is
          incompatible with the modulator type.

        * 'The opt vector must be of length xxx' if opt is used to explicitly
          pass the NTF zeros and these are in the wrong number.

    **Warns:**

        * 'Creating a lowpass ntf.' if the center frequency is different
          from zero, but so low that a low pass modulator must be designed.

        * 'Unable to achieve specified H_inf ...' if the desired H_inf
          cannot be achieved.

        * 'Iteration limit exceeded' if the routine converges too slowly.

    **Notes:**

    This is actually a wrapper function which calls the appropriate version
    of synthesizeNTF, based on the module control flag `optimize_NTF` which
    determines whether to use optimization tools.

    Parameter ``H_inf`` is used to enforce the Lee stability criterion.

    **Example:**

    Fift-order lowpass modulator; zeros optimized for an oversampling ratio of 32.::

        from deltasigma import *
        H = synthesizeNTF(5, 32, 1)
        pretty_lti(H)

    Returns::

              (z -1) (z^2 -1.997z +1) (z^2 -1.992z +0.9999)      
        --------------------------------------------------------
         (z -0.7778) (z^2 -1.796z +0.8549) (z^2 -1.613z +0.665) 

    .. plot::

        from deltasigma import *
        H = synthesizeNTF(5, 32, 1)
        DocumentNTF(H, 32)

    .. seealso::

       :func:`clans` : Closed-Loop Analysis of Noise-Shaper. 
               An alternative method for selecting NTFs based on the 1-norm of the 
               impulse response of the NTF

       :func:`synthesizeChebyshevNTF()` : Select a type-2 highpass Chebyshev NTF.
           This function does a better job than synthesizeNTF if osr
           or H_inf is low.

    """
    if f0 > 0.5:
        raise ValueError('Error. f0 must be less than 0.5.')
    if f0 != 0 and f0 < 0.25/osr:
        warn('Creating a lowpass ntf.')
        f0 = 0
    if f0 != 0 and order % 2 != 0:
        raise ValueError('Order must be even for a bandpass modulator.')
    opt = np.asarray(opt)
    if opt.ndim > 1 or (opt.ndim == 1 and opt.size != order):
        raise ValueError('The opt vector must be of length %d.' % order)

    if optimize_NTF == False:
        ntf = synthesizeNTF0(order, osr, opt, H_inf, f0)
    else:
        ntf = synthesizeNTF1(order, osr, opt, H_inf, f0)
    return ntf

def test_synthesizeNTF():
    """Test function for synthesizeNTF()"""
    from ._utils import cplxpair
    from ._config import optimize_NTF
    optimize_NTF = True
    # synthesizeNTF should have as default values:
    # order=3, osr=64, opt=0, H_inf=1.5, f0=0.0
    z, p, k = synthesizeNTF()
    zref = [1., 1., 1.]
    pref = [.6694, .7654 + .2793j, .7654 - .2793j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # Up next: even order bandpass test
    z, p, k = synthesizeNTF(order=4, osr=32, opt=0, H_inf=1.3, f0=.33)
    zref = [-0.4818 + 0.8763j, -0.4818 - 0.8763j, -0.4818 + 0.8763j,
            -0.4818 - 0.8763j]
    pref = [-0.5125 - 0.7018j, -0.5125 + 0.7018j, -0.3233 - 0.8240j, 
            -0.3233 + 0.8240j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # repeat with zeros optimization
    z, p, k = synthesizeNTF(opt=1)
    zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
    pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # Up next: even order bandpass test
    z, p, k = synthesizeNTF(order=4, osr=32, opt=1, H_inf=1.3, f0=.33)
    zref = [-0.4567 + 0.8896j, -0.4567 - 0.8896j, -0.5064 + 0.8623j,
            -0.5064 - 0.8623j]
    pref = [-0.5125 - 0.7014j, -0.5125 + 0.7014j, -0.3230 - 0.8239j,
            -0.3230 + 0.8239j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    z, p, k = synthesizeNTF(opt=2)
    zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
    pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # Up next: even order bandpass test
    z, p, k = synthesizeNTF(order=4, osr=32, opt=2, H_inf=1.3, f0=.33)
    zref = [-0.4818 + 0.8763j, -0.4818 - 0.8763j, -0.4818 + 0.8763j,
            -0.4818 - 0.8763j]
    pref = [-0.5125 - 0.7018j, -0.5125 + 0.7018j, -0.3233 - 0.8240j,
            -0.3233 + 0.8240j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # opt = 3
    z, p, k = synthesizeNTF(opt=3)
    zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
    pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-3)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-3)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # Up next: even order bandpass test
    z, p, k = synthesizeNTF(order=4, osr=32, opt=3, H_inf=1.3, f0=.33)
    zref = [-0.4567 + 0.8896j, -0.4567 - 0.8896j, -0.5064 + 0.8623j,
            -0.5064 - 0.8623j]
    pref = [-0.5125 - 0.7014j, -0.5125 + 0.7014j, -0.3230 - 0.8239j,
            -0.3230 + 0.8239j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-3)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-3)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # opt = 4
    z, p, k = synthesizeNTF(opt=4)
    zref = [1.0000 + 0.0000j, 0.9993 + 0.0380j, 0.9993 - 0.0380j]
    pref = [0.7652 - 0.2795j, 0.7652 + 0.2795j, 0.6692 + 0.0000j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-3)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-3)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # Up next: odd order lowpass w zero at center band test
    z, p, k = synthesizeNTF(order=5, osr=32, opt=4, H_inf=1.3, f0=0.0)
    zref = [1.0000 + 0.0000j, 0.9986 + 0.0531j, 0.9986 - 0.0531j,
            0.9960 + 0.0892j, 0.9960 - 0.0892j]
    pref = [0.8718 - 0.0840j, 0.8718 + 0.0840j, 0.9390 - 0.1475j,
            0.9390 + 0.1475j, 0.8491 + 0.0000j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-3)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-3)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    # zeros passed explicitly
    opt = [1.0000 + 0.0000j, 0.9986 + 0.06j, 0.9986 - 0.06j,
            0.9960 + 0.0892j, 0.9960 - 0.0892j]
    z, p, k = synthesizeNTF(order=5, osr=32, opt=opt, H_inf=1.3, f0=0.0)
    zref = [1.0000 + 0.0000j, 0.9986 + 0.06j, 0.9986 - 0.06j,
            0.9960 + 0.0892j, 0.9960 - 0.0892j]
    pref = [0.8718 - 0.0840j, 0.8718 + 0.0840j, 0.9390 - 0.1475j,
            0.9390 + 0.1475j, 0.8491 + 0.0000j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-3)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-3)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
