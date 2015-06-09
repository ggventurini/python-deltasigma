# -*- coding: utf-8 -*-
# _simulateQSNR.py
# Module providing the simulateQSNR function
# Copyright 2015 Giuseppe Venturini
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

"""Module providing the simulateQSNR() function
"""

from __future__ import division, print_function

from warnings import warn

import numpy as np
from numpy.fft import fft, fftshift

from ._calculateSNR import calculateSNR
from ._simulateQDSM import simulateQDSM

def simulateQSNR(ntf,
                 R=64,
                 amp=None,
                 f0=0,
                 nlev=2,
                 f=None,
                 k=13):
    """
    Determine the SNR for a quadrature delta-sigma modulator using simulations.

    The modulator is described by a Noise Transfer Function (NTF)
    and the number of quantizer levels.

    Using sine waves located in FFT bins, the SNR is calculated as the ratio
    of the sine wave power to the power in all in-band bins other than those
    associated with the input tone. Due to spectral smearing, the input tone
    is not allowed to lie in bins 0 or 1.

    **Parameters:**

    ntf : tuple, ndarray or LTI object
        The Noise Transfer Function in any form supported by
        :func:`simulateQDSM`, such as an ABCD matrix or an NTF description,
        for example a zpk tuple, num-den tuple or an LTI object.
        If no information is available regarding the STF, it is assumed to
        be unitary.
    R : int, optional
        The oversampling ratio, defining the band of interest. Defaults to 64.
    amp : sequence, optional
        The sequence of the amplitudes to be used for the input signal, in
        ascending order, expressed in dB, where 0 dB means a full-scale (peak
        value :math:`n_{lev}-1`) sine wave. Defaults to [-120 -110...-20 -15 -10
        -9 -8 ... 0] dB.
    f0 : float, optional
        The normalized center frequency of the modulator. Defaults to 0.
    nlev : int, optional
        The number of levels in the modulator quantizer. Defaults to 2.
    f : float, optional
        The input signal frequency, normalized such that :math:`1 \\rightarrow
        f_s`. It is rounded to an FFT bin. If not set, defaults to
        :math:`1/(4\\cdot OSR)`.
    k : int, optional
        The integer ``k`` sets the length of the FFT, which is :math:`2^k`.
        Defaults to 13.

    **Returns:**

    snr : ndarray
        The calculated SNR.
    amp : ndarray
        The amplitude vector corresponding to the SNR.

    """
    if amp is None:
        amp = np.concatenate((np.arange(-120, -20+1, 10),
                              np.atleast_1d(-15),
                              np.arange(-10, 1)))
        if f is None:
            f = f0 + 1./(4*R)
    if np.abs(f - f0) > 1./(2 * R):
        warn('The input tone is out-of-band.')
    N = 2**k
    if N < 8*R:
        # Require at least 8 bins to be "in-band"
        warn('Increasing k to accommodate a large oversampling ratio.')
        k = np.ceil(np.log2(8*R))
        N = 2**k
    F = np.round(f*N)
    if F <= 1:
        warn('Increasing k to accommodate a low input frequency.')
        # We want f*N > 1
        k = np.ceil(np.log2(1./f))
        N = 2**k
        F = 2

    Ntransient = 100
    tone = (nlev - 1)*np.exp(2j*np.pi*F/N*np.arange(-Ntransient, N))
    # Hanning window of length N
    window = 0.5*(1 - np.cos(2*np.pi*np.arange(N)/N))
    f1 = max((np.round(N*(0.5 + f0 - 0.5/R)), 0))
    inBandBins = np.arange(f1, np.round(N*(0.5 + f0 + 0.5/R)) + 1, dtype=np.int32)
    F = F - f1 + N/2.

    snr = np.zeros(amp.shape)
    i = 0
    for A in 10.0**(amp/20.):
        v, _, _, _ = simulateQDSM(A*tone, ntf, nlev)
        hwfft = fftshift(fft(window*v[Ntransient:N + Ntransient]))
        snr[i] = calculateSNR(hwfft[inBandBins], F)
        i = i + 1
    return snr, amp
