# -*- coding: utf-8 -*-
# _synthesizeQNTF.py
# Module providing the synthesizeQNTF function
# Copyright 2013 Giuseppe Venturini
# This file is distributed with python-deltasigma.
#
# python-deltasigma is a 1:1 Python port of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.
#

"""
Module providing the Quadrature NTF synthesis function.
"""

from __future__ import division, print_function

from warnings import warn

import numpy as np

from scipy.signal import cheby2
from scipy.linalg import norm

from ._evalTF import evalTF
from ._dbv import dbv
from ._rmsGain import rmsGain
from ._synthesizeNTF import synthesizeNTF
from ._utils import _get_zpk

# Maximum number of iterations to reach the obj
ITN_MAX = 20

def synthesizeQNTF(order=4, OSR=64, f0=0., NG=-60, ING=-20, n_im=None):
    """Synthesize a noise transfer function for a quadrature modulator.

    **Parameters:**

    order : int, optional
        The order of the modulator. Defaults to 4.

    OSR : int, optional
        The oversampling ratio. Defaults to 64.

    f0 : float, optional
        The center frequency, normalized such that :math:`1 \\rightarrow f_s`.
        Defaults to 0, ie to a low-pass modulator.

    NG : float, optional
        The in-band noise gain, expressed in dB. Defaults to -60.

    ING : float, optional
        The image-band noise gain, in dB. Defaults to -20.

    n_im : int, optional
        The number of in-band image zeros, defaults to ``floor(order/3)``.

    **Returns:**

    ntf : (z, p, k) tuple
        ``ntf`` is a zpk tuple containing the zeros, poles and the gain of the
        synthesized NTF.

    .. note::

        From the Matlab Delta-Sigma Toobox:
        ALPHA VERSION:
        This function uses an experimental ad-hoc method that is
        neither optimal nor robust.

    **Example:**

    Fourth order quadrature modulator::

        from deltasigma import *
        order = 4
        osr = 32
        NG = -50
        ING = -10
        f0 = 1 / 16
        ntf0 = synthesizeQNTF(order, osr, f0, NG, ING)
        pretty_lti(ntf0)

    Returns::

          (z - 0.888 - 0.4598j) (z - 0.9239 + 0.3827j) (z - 0.9239 - 0.3827j) (z - 0.953 - 0.3028j)  
          ---------------------------------------------------------------------------------------------
           (z - 0.5739 - 0.5699j) (z - 0.5913 - 0.2449j) (z - 0.6731 + 0.2788j) (z - 0.8088 - 0.0028j) 


    .. image:: _static/synthesizeQNTF.png


    """
    if n_im is None:
        n_im = np.floor(order/3)
    debug_it = 0
    if n_im == 0:
        # Use synthesizeNTF to get an NTF with the specified NG; ignore ING
        f1 = 0.5/OSR
        x = 1.5
        lowest_f = np.inf
        dfdx = None
        for itn in range(ITN_MAX):
            ntf = synthesizeNTF(order, OSR, 1., x)
            f = dbv(rmsGain(ntf, 0., f1)) - NG
            if debug_it:
                print('x=%.2f f=%.2f' % (x, f))
            if abs(f) < 0.01:
                break
            if dfdx is None:
                dx = 0.1*np.sign(f)
                dfdx = 0
            else:
                dfdx = (f - f_old)/dx
                dx_old = dx
                dx = - f/dfdx
                if abs(dx) > max((1, 2*abs(dx_old))):
                    dx = np.sign(dx)*max((1, 2*abs(dx_old)))
                if x + dx <= 1:
                    # Hinf must be at least 1
                    dx = dx/2.
            f_old = f
            x = x + dx
        if itn == ITN_MAX - 1:
            warn('Iteration limit reached. NTF may be poor.')
        # Rotate the NTF
        z0 = np.exp(2j*np.pi*f0)
        zeros, poles, k = _get_zpk(ntf)
        ntf = (z0*zeros, z0*poles, k)
    else:
        n_in = order - n_im
        f1 = f0 - 0.5/OSR
        f2 = f0 + 0.5/OSR
        z0 = np.exp(2j*np.pi*f0)
        x = np.array([20., 20.])
        # "R" parameters for cheby2()
        lowest_f = np.inf
        dfdx = np.array([float('NaN'), float('NaN')])
        freq = np.linspace(-0.5, 0.5, 200)
        for itn in range(ITN_MAX):
            if debug_it:
                print('\nx = [%.2f, %.2f]' % (x[0], x[1]))
            b1, a1 = cheby2(n_in, x[0], 1./OSR, 'highpass')
            b2, a2 = cheby2(n_im, x[1], 1./OSR, 'highpass')
            ntf0 = (np.concatenate((np.roots(b1)*z0, np.roots(b2)*np.conj(z0))),
                    np.concatenate((np.roots(a1)*z0, np.roots(a2)*np.conj(z0))),
                    1)
            m = evalTF(ntf0, np.exp(2j*np.pi*freq))
            NG0 = dbv(rmsGain(ntf0, f1, f2))
            ING0 = dbv(rmsGain(ntf0, -f1, -f2))
            if debug_it:
                import pylab as plt
                from ._plotPZ import plotPZ
                from ._figureMagic import figureMagic
                plt.figure()
                plt.subplot(121)
                plotPZ(ntf0)
                plt.subplot(122)
                print('NG = %.1f, ING= %.1f' % (NG0, ING0))
                plt.plot(freq, dbv(m))
                figureMagic([-0.5, 0.5], 0.05, 2, [-100, 30], 10, 2)
                plt.hold(True)
                plt.plot([f1, f2], [NG0, NG0], 'k')
                plt.text(np.mean([f1, f2]), NG0, ('NG=%.1fdB' % NG0), va='bottom')
                plt.plot([-f1, -f2], [ING0, ING0], 'k')
                plt.text(np.mean([-f1, -f2]), ING0, ('ING=%.1fdB' % ING0), va='bottom')
                plt.show()
            f = np.array([NG0 - NG, ING0 - ING])
            if max(abs(f)) < 0.01:
                break
            if norm(f) < lowest_f:
                lowest_f = norm(f)
                # ntf0 is ALREADY a zpk tuple
                zeros, poles, k = ntf0
                best = (zeros.copy(), poles.copy(), k)
            if abs(f[0]) > abs(f[1]):
                # adjust x(1)
                i = 0
            else:
                # adjust x(2)
                i = 1
            if np.isnan(dfdx[i]):
                dx = np.sign(f[i])
                dfdx[i] = 0
                dfdx[1 - i] = float('NaN')
            else:
                dfdx[i] = (f[i] - f_old[i])/dx
                dfdx[1 - i] = float('NaN')
                dx = -f[i]/dfdx[i]
                xnew = x[i] + dx
                if xnew < 0.5*x[i]:
                    dx = -0.5*x[i]
                else:
                    if xnew > 2*x[i]:
                        dx = x[i]
            f_old = f.copy()
            x[i] = x[i] + dx
        if itn == ITN_MAX - 1:
            warn('Iteration limit reached. NTF may be poor')
        ntf = best
    return ntf

