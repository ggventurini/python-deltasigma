# -*- coding: utf-8 -*-
# _calculateSNR.py
# Module providing the calculateSNR function
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

"""Module providing the calculateSNR() function
"""

from __future__ import division, print_function
import numpy as np
from numpy.linalg import norm

from ._dbv import dbv

def calculateSNR(hwfft, f, nsig=1):
    """Estimate the SNR from the FFT.

    Estimate the Signal-to-Noise Ratio (SNR), given the in-band bins of
    a Hann-windowed FFT and the location ``f0`` of the input signal (f>0).
    For ``nsig = 1``, the input tone is contained in ``hwfft(f:f+2)``,
    this range is appropriate for a Hann-windowed FFT.

    Each increment in ``nsig`` adds a bin to either side.

    The SNR is expressed in dB.

    **Parameters:**

    hwfft : sequence
            the FFT

    f : integer
        Location of the input signal. Normalized.

    .. note:: f = 0 corresponds to DC, as Python indexing starts from 0.

    nsig : integer, optional
           Extra bins added to either side. Defaults to 1.

    **Returns:**

    SNR : scalar
          The computed SNR value in dB.

    """
    hwfft = hwfft.squeeze()
    signalBins = np.arange(f - nsig + 1, f + nsig + 2, dtype='int64')
    signalBins = signalBins[signalBins > 0]
    signalBins = signalBins[signalBins <= max(hwfft.shape)]
    s = norm(hwfft[signalBins - 1]) # *4/(N*sqrt(3)) for true rms value;
    noiseBins = np.arange(1, max(hwfft.shape) + 1, dtype='int64')
    noiseBins = np.delete(noiseBins, noiseBins[signalBins - 1] - 1)
    n = norm(hwfft[noiseBins - 1])
    if n == 0:
        snr = np.Inf
    else:
        snr = dbv(s/n)
    return snr

