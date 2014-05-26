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

    func:`synthesizeNTF` assumes that magnitude of the denominator of the NTF
    is approximately constant in the passband. When the OSR or ``H_inf`` are
    low, this assumption breaks down and synthesizeNTF yields a non-optimal
    NTF. :func:`synthesizeChebyshevNTF` creates non-optimal NTFs, but fares
    better than synthesizeNTF in the aforementioned circumstances.

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

    **Example:**

    Compare the NTFs created by :func:`synthesizeNTF` and
    :func:`synthesizeChebyshevNTF` when ``OSR`` is low::

        OSR = 4
        order = 8
        H_inf = 3
        H0 = synthesizeNTF(order,OSR,1,H_inf)
        H1 = synthesizeChebyshevNTF(order,OSR,0,H_inf)

    .. plot::

        import pylab as plt
        import numpy as np
        from deltasigma import *
        OSR = 4
        order = 8
        H_inf = 3
        H0 = synthesizeNTF(order,OSR,1,H_inf)
        H1 = synthesizeChebyshevNTF(order,OSR,0,H_inf)
        # 1. Plot the singularities.
        plotsize = (14, 7)
        plt.subplot(121)
        # we plot the singularities of the optimized NTF in light 
        # green with slightly bigger markers so that we can better
        # distinguish the two NTF's when overlayed.
        plotPZ(H1, markersize=7, color='#90EE90')
        plt.hold(True)
        plotPZ(H0, markersize=5)
        plt.title('NTF Poles and Zeros')
        f = np.concatenate((np.linspace(0, 0.75/OSR, 100), np.linspace(0.75/OSR, 0.5, 100)))
        z = np.exp(2j*np.pi*f)
        magH0 = dbv(evalTF(H0, z))
        magH1 = dbv(evalTF(H1, z))
        # 2. Plot the magnitude responses.
        plt.subplot(222)
        plt.plot(f, magH0, label='synthesizeNTF')
        plt.hold(True)
        plt.plot(f, magH1, label='synthesizeChebyshevNTF')
        figureMagic([0, 0.5], 0.05, None, [-80, 20], 10, None, plotsize)
        plt.xlabel('Normalized frequency ($1\\\\rightarrow f_s)$')
        plt.ylabel('dB')
        plt.legend(loc=4)
        plt.title('NTF Magnitude Response')
        # 3. Plot the magnitude responses in the signal band.
        plt.subplot(224)
        fstart = 0.01
        f = np.linspace(fstart, 1.2, 200)/(2*OSR)
        z = np.exp(2j*np.pi*f)
        magH0 = dbv(evalTF(H0, z))
        magH1 = dbv(evalTF(H1, z))
        plt.semilogx(f*2*OSR, magH0, label='synthesizeNTF')
        plt.hold(True)
        plt.semilogx(f*2*OSR, magH1, label='synthesizeChebyshevNTF')
        plt.axis([fstart, 1, -50, 0])
        plt.grid(True)
        sigma_H0 = dbv(rmsGain(H0, 0, 0.5/OSR))
        sigma_H1 = dbv(rmsGain(H1, 0, 0.5/OSR))
        plt.semilogx([fstart, 1], sigma_H0*np.array([1, 1]), linewidth=3, color='#191970')
        plt.text(0.15, sigma_H0 + 1.5, 'RMS gain = %5.0fdB' % sigma_H0)
        plt.semilogx([fstart, 1], sigma_H1*np.array([1, 1]), linewidth=3, color='#228B22')
        plt.text(0.15, sigma_H1 + 1.5, 'RMS gain = %5.0fdB' % sigma_H1)
        plt.xlabel('Normalized frequency ($1\\\\rightarrow f_B$)')
        plt.ylabel('dB')
        plt.legend(loc=3)
        plt.tight_layout()

    Repeat for ``H_inf`` low::

        OSR = 32
        order = 5
        H_inf = 1.2
        H0 = synthesizeNTF(order, OSR, 1, H_inf)
        H1 = synthesizeChebyshevNTF(order, OSR, 1, H_inf)

    .. plot::

        import pylab as plt
        import numpy as np
        from deltasigma import *
        OSR = 32
        order = 5
        H_inf = 1.2
        H0 = synthesizeNTF(order, OSR, 1, H_inf)
        H1 = synthesizeChebyshevNTF(order, OSR, 1, H_inf)
        # 1. Plot the singularities.
        plotsize = (14, 7)
        plt.subplot(121)
        # we plot the singularities of the optimized NTF in light 
        # green with slightly bigger markers so that we can better
        # distinguish the two NTF's when overlayed.
        plotPZ(H1, markersize=7, color='#90EE90')
        plt.hold(True)
        plotPZ(H0, markersize=5)
        plt.title('NTF Poles and Zeros')
        f = np.concatenate((np.linspace(0, 0.75/OSR, 100), np.linspace(0.75/OSR, 0.5, 100)))
        z = np.exp(2j*np.pi*f)
        magH0 = dbv(evalTF(H0, z))
        magH1 = dbv(evalTF(H1, z))
        # 2. Plot the magnitude responses.
        plt.subplot(222)
        plt.plot(f, magH0, label='synthesizeNTF')
        plt.hold(True)
        plt.plot(f, magH1, label='synthesizeChebyshevNTF')
        figureMagic([0, 0.5], 0.05, None, [-80, 20], 10, None, plotsize)
        plt.xlabel('Normalized frequency ($1\\\\rightarrow f_s)$')
        plt.ylabel('dB')
        plt.legend(loc=4)
        plt.title('NTF Magnitude Response')
        # 3. Plot the magnitude responses in the signal band.
        plt.subplot(224)
        fstart = 0.01
        f = np.linspace(fstart, 1.2, 200)/(2*OSR)
        z = np.exp(2j*np.pi*f)
        magH0 = dbv(evalTF(H0, z))
        magH1 = dbv(evalTF(H1, z))
        plt.semilogx(f*2*OSR, magH0, label='synthesizeNTF')
        plt.hold(True)
        plt.semilogx(f*2*OSR, magH1, label='synthesizeChebyshevNTF')
        plt.axis([fstart, 1, -60, -20])
        plt.grid(True)
        sigma_H0 = dbv(rmsGain(H0, 0, 0.5/OSR))
        sigma_H1 = dbv(rmsGain(H1, 0, 0.5/OSR))
        plt.semilogx([fstart, 1], sigma_H0*np.array([1, 1]), linewidth=3, color='#191970')
        plt.text(0.15, sigma_H0 + 1.5, 'RMS gain = %5.0fdB' % sigma_H0)
        plt.semilogx([fstart, 1], sigma_H1*np.array([1, 1]), linewidth=3, color='#228B22')
        plt.text(0.15, sigma_H1 + 1.5, 'RMS gain = %5.0fdB' % sigma_H1)
        plt.xlabel('Normalized frequency ($1\\\\rightarrow f_B$)')
        plt.ylabel('dB')
        plt.legend(loc=3)
        plt.tight_layout()

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
    f_p = None # will be redefined later
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
    from warnings import catch_warnings
    z, p, k = synthesizeChebyshevNTF()
    zref = [1., .9991 + 0.0425j, .9991 - 0.0425j]
    pref = [.6609, .7686 + .2858j, .7686 - .2858j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
    with catch_warnings(record=True) as w:
        z, p, k = synthesizeChebyshevNTF(order=4, OSR=32, opt=1, H_inf=1.5, f0=.33)
        assert len(w) > 0
    zref = [-.4513 + .8924j, -.4513 - .8924j, -.5122 + 0.8589j, -.5122 - 0.8589j]
    pref = [-.2249 + .7665j, -.2249 - .7665j, -.5506 + .6314j, -.5506 - .6314j]
    kref = 1.
    assert np.allclose(cplxpair(z), cplxpair(zref), atol=1e-4, rtol=1e-4)
    assert np.allclose(cplxpair(p), cplxpair(pref), atol=1e-4, rtol=1e-4)
    assert np.allclose(k, kref, atol=1e-4, rtol=1e-4)
