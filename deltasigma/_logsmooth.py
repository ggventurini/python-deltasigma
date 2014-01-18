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
    f1 = int((inBin - 1) % nbin)
    startbin = np.concatenate((np.arange(f1, inBin, nbin), 
                               np.arange(inBin, inBin + 3)
                              ))
    for i in range(1, n + 1):
        startbin = np.concatenate((startbin, 
                       np.arange(startbin[-1] + 1, (i + 1)*inBin, nbin), 
                       (i + 1)*inBin + np.arange(0, 3)
                   ))
    m = startbin[-1] + nbin
    while m < N2:
        startbin = np.concatenate((startbin, np.array((m,))))
        nbin = np.min((nbin*1.1, 2**10))
        m = int(np.round(m + nbin, 0))

    stopbin = np.concatenate((startbin[1:], np.array((N2,)) + 1))
    f = ((startbin + stopbin)/2 - 1)/N
    p = np.zeros(f.shape)
    for i in range(f.shape[0]):
        p[i] = dbp(norm(X[startbin[i]:stopbin[i]])**2/(stopbin[i] - startbin[i]))
    return f, p

def test_logsmooth():
    """Test function for logsmooth()"""
    from ._ds_hann import ds_hann
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
    inBin = np.round(freq/Fs*N)
    fS, US = logsmooth(U, inBin)
    fS = fS[:170]
    US = US[:170]
    fSref, USref = \
     (np.array([  2.16450216e-05,   3.89610390e-05,   5.62770563e-05,
     7.35930736e-05,   9.09090909e-05,   1.08225108e-04,
     1.25541126e-04,   1.42857143e-04,   1.60173160e-04,
     1.77489177e-04,   1.94805195e-04,   2.12121212e-04,
     2.29437229e-04,   2.46753247e-04,   2.64069264e-04,
     2.81385281e-04,   2.98701299e-04,   3.16017316e-04,
     3.33333333e-04,   3.50649351e-04,   3.67965368e-04,
     3.85281385e-04,   4.02597403e-04,   4.19913420e-04,
     4.29653680e-04,   4.31818182e-04,   4.33982684e-04,
     4.36147186e-04,   4.45887446e-04,   4.63203463e-04,
     4.80519481e-04,   4.97835498e-04,   5.15151515e-04,
     5.32467532e-04,   5.49783550e-04,   5.67099567e-04,
     5.84415584e-04,   6.01731602e-04,   6.19047619e-04,
     6.36363636e-04,   6.53679654e-04,   6.70995671e-04,
     6.88311688e-04,   7.05627706e-04,   7.22943723e-04,
     7.40259740e-04,   7.57575758e-04,   7.74891775e-04,
     7.92207792e-04,   8.09523810e-04,   8.26839827e-04,
     8.44155844e-04,   8.58225108e-04,   8.64718615e-04,
     8.66883117e-04,   8.69047619e-04,   8.78787879e-04,
     8.96103896e-04,   9.13419913e-04,   9.30735931e-04,
     9.48051948e-04,   9.65367965e-04,   9.82683983e-04,
     1.00000000e-03,   1.01731602e-03,   1.03463203e-03,
     1.05194805e-03,   1.06926407e-03,   1.08658009e-03,
     1.10389610e-03,   1.12121212e-03,   1.13852814e-03,
     1.15584416e-03,   1.17316017e-03,   1.19047619e-03,
     1.20779221e-03,   1.22510823e-03,   1.24242424e-03,
     1.25974026e-03,   1.27705628e-03,   1.29112554e-03,
     1.29761905e-03,   1.29978355e-03,   1.30194805e-03,
     1.31168831e-03,   1.32900433e-03,   1.34632035e-03,
     1.36363636e-03,   1.38095238e-03,   1.39826840e-03,
     1.41558442e-03,   1.43290043e-03,   1.45021645e-03,
     1.46753247e-03,   1.48484848e-03,   1.50216450e-03,
     1.51948052e-03,   1.53679654e-03,   1.55411255e-03,
     1.57142857e-03,   1.58874459e-03,   1.60606061e-03,
     1.62337662e-03,   1.64069264e-03,   1.65800866e-03,
     1.67532468e-03,   1.69264069e-03,   1.70995671e-03,
     1.72402597e-03,   1.73051948e-03,   1.73268398e-03,
     1.74242424e-03,   1.76082251e-03,   1.78138528e-03,
     1.80411255e-03,   1.82900433e-03,   1.85606061e-03,
     1.88528139e-03,   1.91774892e-03,   1.95346320e-03,
     1.99242424e-03,   2.03571429e-03,   2.08333333e-03,
     2.13528139e-03,   2.19264069e-03,   2.25541126e-03,
     2.32359307e-03,   2.39935065e-03,   2.48268398e-03,
     2.57359307e-03,   2.67424242e-03,   2.78571429e-03,
     2.90800866e-03,   3.04220779e-03,   3.19047619e-03,
     3.35389610e-03,   3.53354978e-03,   3.73051948e-03,
     3.94696970e-03,   4.18506494e-03,   4.44696970e-03,
     4.73593074e-03,   5.05411255e-03,   5.40367965e-03,
     5.78787879e-03,   6.20995671e-03,   6.67424242e-03,
     7.18506494e-03,   7.74675325e-03,   8.36471861e-03,
     9.04437229e-03,   9.79220779e-03,   1.06147186e-02,
     1.15194805e-02,   1.25151515e-02,   1.36103896e-02,
     1.48149351e-02,   1.61396104e-02,   1.75974026e-02,
     1.92012987e-02,   2.09653680e-02,   2.29058442e-02,
     2.50303030e-02,   2.72467532e-02,   2.94632035e-02,
     3.16796537e-02,   3.38961039e-02,   3.61125541e-02,
     3.83290043e-02,   4.05454545e-02]), np.array([ 
    -3.38256602e+02,  -3.35680802e+02,  -3.34955637e+02,
    -3.31327127e+02,  -3.25855874e+02,  -3.29103355e+02,
    -3.20949679e+02,  -3.32450849e+02,  -3.35225202e+02,
    -3.34485568e+02,  -3.30257414e+02,  -3.29290454e+02,
    -3.19938944e+02,  -3.26303196e+02,  -3.18438211e+02,
    -3.08049597e+02,  -3.06579645e+02,  -3.18589438e+02,
    -3.22795585e+02,  -3.20720313e+02,  -3.23943465e+02,
    -3.18410178e+02,  -3.14006703e+02,  -3.06887520e+02,
    -6.01408794e+00,   6.51197554e-03,  -6.01408794e+00,
    -2.97355915e+02,  -3.20713161e+02,  -3.14764659e+02,
    -3.19601722e+02,  -3.25738958e+02,  -3.24901490e+02,
    -3.29780263e+02,  -3.30378006e+02,  -3.26494435e+02,
    -3.30794707e+02,  -3.20546818e+02,  -3.30534642e+02,
    -3.19514346e+02,  -3.34461369e+02,  -3.32145601e+02,
    -3.29824120e+02,  -3.30165209e+02,  -3.31177740e+02,
    -3.22895278e+02,  -3.34922380e+02,  -3.27925810e+02,
    -3.30624932e+02,  -3.33688205e+02,  -3.33118130e+02,
    -3.21601327e+02,  -5.55090747e+01,  -4.24987747e+01,
    -4.85193746e+01,  -3.34050272e+02,  -3.29761613e+02,
    -3.32017540e+02,  -3.27930566e+02,  -3.30543187e+02,
    -3.12454566e+02,  -3.27485700e+02,  -3.20772849e+02,
    -3.17857749e+02,  -3.14837482e+02,  -3.19784038e+02,
    -3.26885855e+02,  -3.30049777e+02,  -3.26864139e+02,
    -3.25721205e+02,  -3.19933092e+02,  -3.16392632e+02,
    -3.05698899e+02,  -3.11597018e+02,  -3.20731465e+02,
    -3.23881857e+02,  -3.24672373e+02,  -3.29662306e+02,
    -3.28122258e+02,  -3.32338288e+02,  -8.50514998e+01,
    -7.20411998e+01,  -7.80617997e+01,  -3.39815014e+02,
    -3.34795713e+02,  -3.33788288e+02,  -3.37592846e+02,
    -3.34135813e+02,  -3.32982286e+02,  -3.34403852e+02,
    -3.38458853e+02,  -3.32145223e+02,  -3.34008178e+02,
    -3.18018946e+02,  -3.36678711e+02,  -3.35114520e+02,
    -3.34636076e+02,  -3.28243681e+02,  -3.24123154e+02,
    -3.07175476e+02,  -3.19669353e+02,  -3.24411950e+02,
    -3.31757229e+02,  -3.27256233e+02,  -3.33953356e+02,
    -3.28750172e+02,  -3.34046831e+02,  -3.22633970e+02,
    -7.70926996e+01,  -6.40823997e+01,  -7.01029996e+01,
    -3.30474657e+02,  -3.32634583e+02,  -3.29125133e+02,
    -3.30037433e+02,  -3.29135335e+02,  -3.23031499e+02,
    -3.15327586e+02,  -3.28033692e+02,  -3.27416910e+02,
    -3.23708823e+02,  -3.34949826e+02,  -3.33223824e+02,
    -3.26295841e+02,  -3.20128696e+02,  -3.33775118e+02,
    -3.24591252e+02,  -3.16785294e+02,  -3.16060771e+02,
    -3.35340249e+02,  -3.36352885e+02,  -3.31545433e+02,
    -3.35918671e+02,  -3.31247086e+02,  -3.22353315e+02,
    -3.30659242e+02,  -3.26197099e+02,  -3.33943174e+02,
    -3.27104794e+02,  -3.28393533e+02,  -3.25787147e+02,
    -3.27725773e+02,  -3.28354321e+02,  -3.33420503e+02,
    -3.29958863e+02,  -3.33004059e+02,  -3.27341144e+02,
    -3.32645042e+02,  -3.26692431e+02,  -3.22999148e+02,
    -3.22628391e+02,  -3.22379886e+02,  -3.00321481e+02,
    -3.00705013e+02,  -3.22941561e+02,  -3.25563520e+02,
    -3.27831432e+02,  -3.25586349e+02,  -3.27610733e+02,
    -3.20556431e+02,  -3.06655655e+02,  -3.07172708e+02,
    -3.21314379e+02,  -3.21324902e+02,  -3.20683384e+02,
    -3.12945454e+02,  -3.13687003e+02,  -3.30837064e+02,
    -3.26718909e+02,  -3.07344944e+02]))
    assert np.allclose(fS, fSref, atol=1e-8, rtol=1e-5)
    assert np.allclose(fS, fSref, atol=1e-8, rtol=1e-5)
