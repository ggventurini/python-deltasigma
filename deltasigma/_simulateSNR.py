# -*- coding: utf-8 -*-
# _simulateSNR.py
# Module providing the simulateSNR function
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

"""Module providing the simulateSNR() function
"""

from __future__ import division
import collections
from warnings import warn
import numpy as np
from numpy.fft import fft, fftshift

from ._simulateDSM import simulateDSM
from ._calculateSNR import calculateSNR
from ._mapQtoR import mapQtoR
from ._utils import _get_zpk

def simulateSNR(arg1, osr, amp=None, f0=0, nlev=2, f=None, k=13,
                quadrature=False):
    """Determine the SNR for a delta-sigma modulator by using simulations.

    Simulate a delta-sigma modulator with sine wave inputs of various
    amplitudes and calculate the signal-to-noise ratio (SNR) in dB for each 
    input.

    Three alternative descriptions of the modulator can be used:

     * The modulator is described by a noise transfer function (NTF), provided
       as ``arg1`` and the number of quantizer levels (``nlev``).

     * Alternatively, the first argument to simulateSNR may be an ABCD matrix. 

     * Lastly, ``arg1`` may be a function taking the input signal as its
       sole argument.

    The band of interest is defined by the oversampling ratio (``osr``)
    and the center frequency (``f0``).

    The input signal is characterized by the ``amp`` vector and the ``f`` variable.
    A default value for ``amp`` is used if not supplied.

    ``f`` is the input frequency, normalized such that 1 -> fs;
    ``f`` is rounded to an FFT bin.

    Using sine waves located in FFT bins, the SNR is calculated as the ratio
    of the sine wave power to the power in all in-band bins other than those
    associated with the input tone. Due to spectral smearing, the input tone
    is not allowed to lie in bins 0 or 1. The length of the FFT is :math:`2^k`.

    If the NTF is complex, :func:`simulateQDSM` (which is slow, also available
    in a future release) is called.

    If ABCD is complex, :func:`simulateDSM` is used with the real equivalent
    of ABCD in order to speed up simulations.

    Future versions may accommodate STFs.

    **Parameters:**

    arg1 : scipy 'lti' object, or ndarray
        The first argument may be one of the various supported representations
        for a (SISO) transfer function or an ABCD matrix.

    osr : int
        The over-sampling ratio.

    amp : sequence, optional
        The amplitudes in dB, referred to the FS, for which the SNR is to be 
        evaluated. ``amp`` defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0]dB,
        where 0 dB means a full-scale (peak value = nlev-1) sine wave.

    f0 : float, optional
        The center frequency. Normalized. Defaults to 0.

    nlev : int, optional
        Number of quantizer levels, defaults to 2.

    f : float, optional
        Test signal input frequency. Normalized. Rounded to an FFT bin.
        Defaults to:

    .. math::

        f = \\frac{1}{4\ \mathrm{OSR}}

    k : int, optional
        The number of samples used to compute the FFT is set by the integer `k`
        - default value 13 - through:

    .. math::

        N_{samples} = 2^k

    quadrature : boolean, optional
        Whether the delta sigma modulator is a quadrature modulator or not.
        Defaults to ``False``.

    .. note:: Setting ``quadrature`` to ``True`` results in a \\
        ``NotImplementedError`` being raised, as :func:`simulateQDSM` has not been \\
        implemented yet.

    **Returns:**

    snr : ndarray
        The SNR, from simulation.

    amp : ndarray
        The amplitudes corresponding to the SNR values.

    .. rubric:: Example:

    Compare the SNR vs input amplitude curve for a fifth-order modulator, as 
    determined by the describing function method (:func:`predictSNR`) with
    that determined by simulation (:func:`simulateSNR`).::

        import pylab as plt
        from deltasigma import *
        OSR = 32
        H = synthesizeNTF(5, OSR, 1)
        snr_pred, amp, _, _, _ = predictSNR(H,OSR)
        snr, amp = simulateSNR(H, OSR)
        plt.plot(amp, snr_pred, 'b', label='Predicted')
        plt.hold(True)
        plt.plot(amp, snr, 'go', label='Simulated')
        plt.grid(True)
        figureMagic([-100, 0], 10, None,
                    [0, 100], 10, None)
        plt.xlabel('Input Level, dB')
        plt.ylabel('SNR, dB')
        s = 'peak SNR = %4.1fdB\\n' % max(snr)
        plt.text(-65, 15, s, horizontalalignment='left')
        plt.legend(loc='best')

    .. plot::

        import pylab as plt
        from deltasigma import *
        OSR = 32
        H = synthesizeNTF(5, OSR, 1)
        snr_pred, amp, _, _, _ = predictSNR(H,OSR)
        snr, amp = simulateSNR(H, OSR)
        plt.plot(amp, snr_pred, 'b', label='Predicted')
        plt.hold(True)
        plt.plot(amp, snr, 'go', label='Simulated')
        plt.grid(True)
        figureMagic([-100, 0], 10, None,
                    [0, 100], 10, None)
        plt.xlabel('Input Level, dB')
        plt.ylabel('SNR, dB')
        s = 'peak SNR = %4.1fdB\\n' % max(snr)
        plt.text(-65, 15, s, horizontalalignment='left')
        plt.legend(loc='best')

    """
    # Look at arg1 and decide if the system is quadrature
    quadrature_ntf = False
    if callable(arg1):
        pass
    elif hasattr(arg1, '__class__') and arg1.__class__.__name__ == 'lti':
        # scipy LTI object
        for roots in _get_zpk(arg1)[:2]:
            if np.any(np.abs(np.imag(np.poly(roots))) > 0.0001):
                quadrature = True
                quadrature_ntf = True
    else: # ABCD matrix
        if not np.all(np.imag(arg1) == 0):
            quadrature = True

    if amp is None:
        amp = np.concatenate((
                              np.arange(- 120, -20 + 1, 10),
                              np.array((-15,)),
                              np.arange(-10, 1)
                            ))
    elif not isinstance(amp, collections.Iterable):
        amp = np.array((amp, ))
    else:
        amp = np.asarray(amp)
    osr_mult = 2
    if f0 != 0 and not quadrature:
        osr_mult = 2*osr_mult
    if f is None or np.isnan(f):
        f = f0 + 0.5/(osr*osr_mult) # Halfway across the band
    M = nlev - 1
    if quadrature and not quadrature_ntf:
        # Modify arg1 (ABCD) and nlev so that simulateDSM can be used
        nlev = np.array([nlev, nlev]).reshape(1, -1)
        arg1 = mapQtoR(arg1)
    if abs(f - f0) > 1/(osr*osr_mult):
        warn('The input tone is out-of-band.')
    N = 2**k
    if N < 8*2*osr:
        warn('Increasing k to accommodate a large oversampling ratio.')
        k = np.array(np.ceil(np.log2(8*2*osr)), dtype=np.int64)
        N = 2**k
    F = np.round(f*N)
    if np.abs(F) <= 1:
        warn('Increasing k to accommodate a low input frequency.')
        # Want f*N >= 1
        k = np.ceil(np.log2(1./f))
        N = 2**k
        F = 2
    Ntransient = 100
    soft_start = 0.5*(1 - np.cos(2*np.pi/Ntransient * \
                                 np.arange(0, Ntransient/2)))
    if not quadrature:
        tone = M*np.sin(2*np.pi*F/N*np.arange(0, N + Ntransient))
        tone[:Ntransient/2] = tone[:Ntransient/2] * soft_start
    else:
        tone = M * np.exp(2j*np.pi*F/N * np.arange(0, N + Ntransient))
        tone[:Ntransient/2] = tone[:Ntransient/2] * soft_start
        if not quadrature_ntf:
            tone = np.hstack((np.real(tone), np.imag(tone)))
    # create a Hann window
    window = 0.5*(1 - np.cos(2*np.pi*np.arange(0, N)/N))
    if quadrature:
        window = np.vstack((window, window))
    if f0 == 0:
        inBandBins = int(N/2) + np.arange(3,
                                     np.round(N/osr_mult/osr) + 1,
                                     dtype=np.int32)
        F = F - 2
    else:
        f1 = np.round(N*(f0 - 1./osr_mult/osr))
        # Should exclude DC
        inBandBins = int(N/2) + np.arange(f1,
                                     np.round(N*(f0 + 1./osr_mult/osr)) + 1,
                                     dtype=np.int32)
        F = F - f1 + 1
    snr = np.zeros(amp.shape)
    i = 0
    for A in np.power(10.0, amp/20):
        if callable(arg1):
            v = arg1(A*tone)
        else:
            if quadrature_ntf:
                raise NotImplementedError("simulateQDSM has not been \
                                           implemented yet.")
                #v = simulateQDSM(A*tone, arg1, nlev)
            else:
                v, _, _, _ = simulateDSM(A*tone, arg1, nlev)
                if quadrature:
                    v = v[0, :] + 1j*v[1, :]
        hwfft = fftshift(fft(window*v[Ntransient:N + Ntransient]))
        snr[i] = calculateSNR(hwfft[inBandBins - 1], F)
        i += 1
    return snr, amp

