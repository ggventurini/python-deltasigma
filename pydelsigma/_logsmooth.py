# -*- coding: utf-8 -*-
# _logsmooth.py
# Module providing the logsmooth function
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

"""Module providing the logsmooth() function
"""

from __future__ import division
import numpy as np
from scipy.linalg import norm

from ._dbp import dbp

def logsmooth(X, inBin, nbin=8, n=3):
    """Smooth the fft, X, and convert it to dB.

    Use nbin bins from 0 to 3*inBin, 
    thereafter increase bin sizes by a factor of 1.1, staying less than 2^10.

    For the n sets of bins inBin+[0:2], 2*inBin+[0:2], ... 
    n*inBin+[0:2], don't do averaging. This way, the noise BW
    and the scaling of the tone and its harmonics are unchanged.

    Unfortunately, harmonics above the nth appear smaller than they 
    really are because their energy is averaged over many bins.
    """
    if not np.prod(X.shape) == max(X.shape):
        raise ValueError('Expected a (N,) or (N, 1)-shaped array.')
    if len(X.shape) > 1:
        X = np.squeeze(X)
    N = X.shape[0]
    N2 = np.floor(N/2)
    f1 = ((inBin - 1) % nbin) + 1
    startbin = np.concatenate((np.arange(f1, inBin, nbin), 
                               np.arange(inBin, inBin + 3)
                              ))
    for i in range(n):
        startbin = np.concatenate((startbin, 
                       np.arange(startbin[-1] + 1, (i + 1)*inBin, nbin), 
                       (i + 1)*inBin + np.arange(0, 3)
                   ))
    m = startbin[-1] + nbin
    while m < N2:
        startbin = np.concatenate((startbin, np.array((m,))))
        nbin = np.min((nbin*1.1, 2**10))
        m = int(np.round(m + nbin, 0))

    stopbin = np.concatenate((startbin[1:] - 1, np.array((N2,))))
    f = ((startbin + stopbin)/2 - 1)/N
    p = np.zeros(f.shape)
    for i in range(max(f.shape) + 1):
        p[i] = dbp(norm(X[startbin[i]:stopbin[i]])**2/(stopbin[i] - startbin[i] + 1))
    return f, p

