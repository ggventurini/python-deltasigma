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
    """Smooth the fft, and convert it to dB.

    **Parameters:**

    X : (N,) ndarray
        The FFT data.

    inBin : int
        The bin index of the input sine wave (if any).

    nbin : int, optional
        The number of bins on which the averaging will be performed,
        used *before* 3*inBin

    n : int, optional
        Around the location of the input signal and its harmonics (up to the
        third harmonic), don't average for n bins.

    The logsmooth algorithm uses nbin bins from 0 to 3*inBin,
    thereafter the bin sizes is increased by a factor of 1.1,
    staying less than 2^10.

    For the :math:`n` sets of bins:
    :math:`inBin + i, 2*inBin + i ... n*inBin+i`, where :math:`i \\in [0,2]`
    don't do averaging. This way, the noise BW
    and the scaling of the tone and its harmonics are unchanged.

    .. note::

        Unfortunately, harmonics above the nth appear smaller than they
        really are because their energy is averaged over many bins.

    **Returns:**

    f, p : tuple of 1d- ndarrays
        The bins and smoothed FFT, expressed in dB.

    .. seealso::

         * :func:`plotSpectrum`, convenience function to first call
           :func:`logsmooth` and then plot on a logarithmic x-axis its return
           value.

         * :func:`circ_smooth`, smoothing algorithm suitable for linear
           x-axis plotting.

    .. plot::

        import matplotlib.pyplot as plt
        import numpy as np
        from deltasigma import dbv, ds_hann, figureMagic, logsmooth
        T = 2 #s
        Fs = 231e3 #Hz
        N = int(np.round(T*Fs, 0)) # FFT points
        freq = .1e3
        t = np.arange(N)/Fs
        u0 = np.sin(2*np.pi*t*freq)
        u0 = u0 + .01*u0**2+.001*u0**3+.005*u0**4
        U = np.fft.fft(u0 * ds_hann(N))/(N/4)
        f = np.linspace(0, Fs, N + 1)
        f = f[:N/2 + 1]
        plt.subplot(211)
        plt.semilogx(f, dbv(U[:N/2 + 1]))
        plt.hold(True)
        inBin = np.round(freq/Fs*N)
        fS, US = logsmooth(U, inBin)
        plt.semilogx(fS*Fs, US, 'r', linewidth=2.5)
        plt.xlim([f[0]*Fs, Fs/2])
        plt.ylabel('U(f) [dB]')
        figureMagic(xRange=[100, 1e4], yRange=[-400, 0], name='Spectrum')
        plt.subplot(212)
        plt.loglog(fS[1:]*Fs, np.diff(fS*Fs))
        plt.xlabel('f [Hz]')
        plt.ylabel('Averaging interval [Hz]')
        figureMagic(xRange=[100, 1e4])
        plt.show()

    """
    # preliminary sanitization of the input
    if not np.prod(X.shape) == max(X.shape):
        raise ValueError('Expected a (N,) or (N, 1)-shaped array.')
    if len(X.shape) > 1:
        X = np.squeeze(X)
    inBin = int(inBin)

    N = X.shape[0]
    N2 = int(np.floor(N/2))
    f1 = int(inBin % nbin)
    startbin = np.concatenate((np.arange(f1, inBin, nbin), 
                               np.arange(inBin, inBin + 3)
                              ))
    i = 1 # my fix
    while i < n: # n can be big and xrange is not in Python3
        startbin = np.concatenate((startbin, 
                       np.arange(startbin[-1] + 1, (inBin + 1)*(i + 1) - 1, nbin), 
                       (i + 1)*(inBin + 1) - 1 + np.arange(0, 3)
                   ))
        i = i + 1
    startbin = np.concatenate((startbin, np.array((startbin[-1] + 1,)))) # my fix
    m = startbin[-1] + nbin
    while m < N2 - 1:
        startbin = np.concatenate((startbin, np.array((m,))))
        nbin = np.min((nbin*1.1, 2**10))
        m = int(np.round(m + nbin, 0))

    stopbin = np.concatenate((startbin[1:] - 1, np.array((N2 - 1,))))
    f = ((startbin + stopbin)/2)/N
    p = np.zeros(f.shape)
    for i in range(f.shape[0]):
        p[i] = dbp(norm(X[startbin[i]:stopbin[i] + 1])**2/(stopbin[i] - startbin[i] + 1))
    return f, p

def test_logsmooth():
    """Test function for logsmooth()"""
    import pkg_resources, scipy.io
    from ._dbv import dbv
    from ._ds_f1f2 import ds_f1f2
    from ._ds_hann import ds_hann
    from ._simulateDSM import simulateDSM
    from ._synthesizeNTF import synthesizeNTF
    from ._undbv import undbv
    from numpy.fft import fft, fftshift
    import numpy as np

    f0 = 1./8
    OSR = 64
    order = 8
    N = 8192
    H = synthesizeNTF(order, OSR, 1, 1.5, f0)
    fB = int(np.ceil(N/(2. * OSR)))
    quadrature = False
    M = 1
    f1, f2 = ds_f1f2(OSR, f0, quadrature)
    delta = 2
    Amp = undbv(-3)
    f1_bin = np.round(f1*N)
    f2_bin = np.round(f2*N)
    f = .3
    fin = np.round(((1 - f)/2*f1 + (f + 1)/2 * f2) * N)
    t = np.arange(0, N).reshape((1, -1))
    u = Amp * M * np.cos((2*np.pi/N)*fin*t)
    v, xn, xmax, y = simulateDSM(u, H, M + 1)

    window = ds_hann(N)
    NBW = 1.5/N
    spec0 = fft(v * window) / (M*N/4)
    freq = np.linspace(0, 0.5, N/2 + 1)

    fl, pl = logsmooth(spec0, fin - 1) ### -1 is really important THIS IS NOT MATLAB!!

    fname = pkg_resources.resource_filename(__name__, "test_data/test_logsmooth.mat")
    data = scipy.io.loadmat(fname)
    assert np.allclose(fl, data['fl'], atol=1e-8, rtol=1e-5)
    assert np.allclose(pl, data['pl'], atol=1e-8, rtol=1e-5)
