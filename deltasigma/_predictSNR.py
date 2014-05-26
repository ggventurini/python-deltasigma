# -*- coding: utf-8 -*-
# _predictSNR.py
# Module providing predictSNR
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

"""Module providing the predictSNR function.
"""


from __future__ import division
import numpy as np
from scipy.signal import dimpulse, freqz
from scipy.special import erfinv
from scipy.interpolate import interp1d

from ._dbp import dbp
from ._utils import _get_num_den

def predictSNR(ntf, R=64, amp=None, f0=0.):
    """Predict the SNR curve of a binary delta-sigma modulator.

    The prediction is performed using the describing function method of Ardalan
    and Paulos [2]_ .

    **Parameters:**
    
    ntf : lti or tf object, or zpk or (num, den) tuples
        The noise transfer function specifying the modulator.
    
    R : scalar, optional
        The oversampling ratio, defaults to 64.
    
    amp : ndarray-like, optional
        The magnitudes to be used for the input signal. They are expressed in 
        dB, where 0 dB means a full-scale (peak value = 1) sine wave. 
        Defaults to [-120 -110...-20 -15 -10 -9 -8 ... 0].
    
    f0 : scalar, 
        The normalized input signal frequency. Defaults to 0.
    
    **Notes:**
    
    The band of interest is defined by the oversampling ratio (``R``) and the 
    center frequency (``f0``).

    The algorithm assumes that the ``amp`` vector is sorted in increasing order;
    once instability is detected, the remaining SNR values are set to ``-Inf``.
    
    Future versions may accommodate STFs.

    **Returns:**

    snr : ndarray
        A vector of SNR values, in dB.

    amp : ndarray
        A vector of amplitudes, in dB.

    k0 : ndarray
        The quantizer signal gain.

    k1: ndarray
        The quantizer noise gain.

    sigma_e2 : scalar
        The power of the quantizer noise (not in dB).
    

    .. rubric:: Implementation details:

    The describing function method of A&P treats the quantizer processes
    signal and noise components separately. The quantizer is modelled as two
    (not necessarily equal) linear gains, :math:`k_0` (``k0`` in the code) 
    and :math:`k_1` (``k1``), and an additive white gaussian noise source of
    power :math:`\\sigma_e^2` (``sigma_e2``), as shown in the figure below. 

    :math:`k_0`, :math:`k_1` and :math:`\\sigma_e^2` are calculated as
    functions of the input.

    .. image:: ../doc/_static/predictSNR.png
        :align: center
        :alt: modulator model for predictSNR


    The modulator's loop filter is assumed to have nearly infinite gain at
    the test frequency.

    .. rubric:: Example:

    See :func:`simulateSNR` for an example use of this function.

    .. rubric:: References

    .. [2] Ardalan, S.H.; Paulos, J.J., "An analysis of nonlinear behavior in
           delta - sigma modulators," Circuits and Systems, IEEE Transactions
           on, vol.34, no.6, pp.593,603, Jun 1987
    
    """

    # extract num, den
    if (hasattr(ntf, 'inputs') and not ntf.inputs == 1) or \
       (hasattr(ntf, 'outputs') and not ntf.outputs == 1):
        raise TypeError("The supplied TF isn't a SISO transfer function.")

    num, den = _get_num_den(ntf)
    Nb = 100
    if f0 == 0:
        band_of_interest = np.linspace(0, np.pi/R, Nb)
    else:
        band_of_interest = np.linspace(2*np.pi*(f0 - 0.25/R), 2*np.pi*(f0 + 0.25/R), Nb)
        XTAB = np.linspace(-2, 0, 21)
        YTAB = np.array([
            [0.46575960516930,   0.67366999387741],
            [0.47904652357101,   0.68426650762558],
            [0.49316295981407,   0.69527947902679],
            [0.50817364454269,   0.70673173666000],
            [0.52414894104004,   0.71864765882492],
            [0.54116523265839,   0.73105299472809],
            [0.55930554866791,   0.74397552013397],
            [0.57866013050079,   0.75744456052780],
            [0.59932720661163,   0.77149158716202],
            [0.62141352891922,   0.78615015745163],
            [0.64503526687622,   0.80145609378815],
            [0.67031890153885,   0.81744754314423],
            [0.69740217924118,   0.83416539430618],
            [0.72643494606018,   0.85165339708328],
            [0.75758063793182,   0.86995816230774],
            [0.79101717472076,   0.88912981748581],
            [0.82693856954575,   0.90922164916992],
            [0.86555624008179,   0.93029111623764],
            [0.90710091590881,   0.95239937305450],
            [0.95182400941849,   0.97561222314835],
            [1.00000000000000,   1.00000000000000]])

    if amp is None:
        amp = np.concatenate((np.arange(- 120, -20 + 1, 10),
                              np.array((-15,)),
                              np.arange(-10, 1)
                            ))
    num = np.real_if_close(num)
    den = np.real_if_close(den)
    num1 = num - den
    N = max(amp.shape)
    snr = np.zeros((1, N)) - np.Inf
    k0 = np.zeros((1, N))
    k1 = np.zeros((1, N))
    sigma_e2 = np.zeros((1, N))
    u = 10.0**(amp/20)
    Nimp = 100
    unstable = False
    for n in range(N):
        # Calculate sigma_e2
        if f0 == 0:
            erfinvu = erfinv(u[n])
            sigma_e2[0, n] = 1 - u[n]**2 - 2/np.pi * np.exp(-2*erfinvu**2)
        else:
            # % Sinusoidal input.
            # Solve sqrt(pi)*u/2 = rho * hypergeo(0.5,2,-rho^2);
            # Formulate as solve f(rho) = 0, f = rho*M(0.5,2,-rho^2)-K
            # and use the secant method.
            K = 0.5*np.sqrt(np.pi)*u[n]
            if n == 0:
                # Initial guess; otherwise use previous value.
                rho = u[n]**2
                fprime = 1
            drho = 1
            f_prev = None
            for itn in range(0, 20):
                m0 = interp1d(XTAB, YTAB[:, 1], kind='cubic')(-rho**2)
                f = rho*m0 - K
                if itn > 0:
                    fprime = max((f - f_prev)/drho, 0.5) #Secant approx.
                if abs(f) < 1e-08:
                    break #!Converged
                drho = -f/fprime
                if abs(drho) > 0.2:
                    drho = np.sign(drho) * 0.2
                if abs(drho) < 1e-06:
                    break #!Converged
                rho = rho + drho
                f_prev = f
            m1 = interp1d(XTAB, YTAB[:, 0], kind='cubic')(-rho**2)
            sigma_e2[0, n] = 1 - u[n]**2/2 - 2/np.pi*m1**2
        # Iterate to solve for k1 and sigma_1.
        # Using one of MATLAB's nonlinear equation solvers would be more efficient,
        # but this function code would then require the optimization toolbox.
        # !Future work: put in 2-D BFGS code.
        if n > 0:
            #  Use the previous value of k1 as the initial guess.
            k1[0, n] = k1[0, n - 1]
        else:
            k1[0, n] = 1.2
        k1_prev = 0
        itn = 0
        if f0 == 0:
            k1sigma1 = np.sqrt(2/np.pi) * np.exp(-erfinvu**2)
        else:
            k1sigma1 = np.sqrt(2/np.pi)*m1
        while abs(k1[0, n] - k1_prev) > 1e-06*(1 + k1[0, n]) and itn < 100:
            #  Create the function: H_hat = L1/(1-k1*L1)=(H-1)/(H*(1-k1)+k1).
            den1 = (1 - k1[0, n])*num + den*k1[0, n]
            #  Calculate pGain, the square of the 2-norm of H_hat.
            pGain, Nimp = powerGain(num1, den1, Nimp)
            if np.isinf(pGain):
                unstable = True
                break
            sigma_1 = np.sqrt(pGain * sigma_e2[0, n])
            k1_prev = k1[0, n]
            k1[0, n] = k1sigma1/sigma_1
            itn = itn + 1

        if unstable:
            break
        if f0 == 0:
            y0 = np.sqrt(2)*erfinvu*sigma_1
            k0[0, n] = u[n]/y0
        else:
            k0[0, n] = np.sqrt(2/np.pi)*m0/sigma_1
        _, h = freqz(num, (1 - k1[0, n])*num + k1[0, n]*den, band_of_interest)
        # For both DC and sine wave inputs, use u^2/2 as the signal 
        # power since true DC measurements are usually impossible.
        snr[0, n] = dbp(0.5*u[n]**2/(np.sum(h**2)/(R*Nb)*sigma_e2[0, n]))
    return snr.squeeze(), amp.squeeze(), k0.squeeze(), k1.squeeze(), sigma_e2.squeeze()

def powerGain(num, den, Nimp=100):
    """Calculate the power gain of a TF given in coefficient form.

    Nimp is the recommended number of impulse response samples for use
    in future calls and Nimp0 is the suggested number (100) to use.
    """
    unstable = False
    _, (imp, ) = dimpulse((num, den, 1), t=np.linspace(0, Nimp, Nimp))
    if np.sum(abs(imp[Nimp - 11:Nimp])) < 1e-08 and Nimp > 50:
        Nimp = np.round(Nimp/1.3)
    else:
        while np.sum(abs(imp[Nimp - 11:Nimp])) > 1e-06:
            Nimp = Nimp*2
            _, (imp, ) = dimpulse((num, den, 1), t=np.linspace(0, Nimp, Nimp))
            if np.sum(abs(imp[Nimp - 11:Nimp])) >= 50 or Nimp >= 10000.0:
                unstable = True
                break

    if not unstable:
        pGain = np.sum(imp**2)
    else:
        pGain = np.Inf

    return pGain, Nimp

def test_predictSNR():
    """Test function for predictSNR()"""
    from ._synthesizeNTF import synthesizeNTF
    snr_ref = [-32.3447, -22.3447,  -12.3447,   -2.3447,    7.6553,   17.6553,
               27.6553,   37.6552,   47.6545,   57.6475,   67.5768,   72.4043,
               76.8266,   77.5913,   78.2773,   78.8451,   79.2116,   79.0974,
               -np.Inf,      -np.Inf,      -np.Inf,      -np.Inf,      -np.Inf]
    amp_ref = [-120,  -110,  -100,   -90,   -80,   -70,   -60,   -50,   -40,   -30,   -20,
               -15,    -10,    -9,    -8,    -7,    -6,    -5,    -4,    -3,    -2,    -1,
               +0]
    k0_ref = [1.6289,    1.6289,    1.6289,    1.6289,    1.6289,    1.6289,
              1.6289,    1.6289,    1.6288,    1.6283,    1.6227,    1.6088,
              1.5596,    1.5383,    1.5088,    1.4663,    1.4003,    1.2747,
              0,         0,         0,         0,         0]

    k1_ref = [1.6289,    1.6289,    1.6289,    1.6289,    1.6289,    1.6289,
              1.6289,    1.6289,    1.6287,    1.6274,    1.6142,    1.5819,
              1.4752,    1.4326,    1.3768,    1.3025,    1.1995,    1.0387,
              0.3706,         0,         0,         0,         0]

    sigma_e2_ref = [0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3631,
                    0.3607,    0.3591,    0.3566,    0.3525,    0.3459,    0.3352,
                    0.3178,         0,         0,         0,         0]

    order = 5
    osr = 32
    f0 = 0
    Hinf = 1.5

    ntf = synthesizeNTF(order, osr, 2, Hinf, f0)
    snr_pred, amp_pred, k0, k1, sigma_e2 = predictSNR(ntf, osr, None, f0)

    # Delete garbage data
    for check in (np.isinf, np.isnan):
        i = check(snr_ref)
        snr_ref = np.delete(snr_ref, np.where(i))
        amp_ref = np.delete(amp_ref, np.where(i))
        i = check(snr_pred)
        snr_pred = np.delete(snr_pred, np.where(i))
        amp_pred = np.delete(amp_pred, np.where(i))

    assert np.allclose(snr_pred, snr_ref, atol=1e-2, rtol=1e-3)
    assert np.allclose(amp_pred, amp_ref, atol=1e-2, rtol=5e-4)
    assert np.allclose(k0, k0_ref, atol=1e-3, rtol=1e-2)
    assert np.allclose(k1, k1_ref, atol=1e-3, rtol=50e-2)
    assert np.allclose(sigma_e2, sigma_e2_ref, atol=1e-3, rtol=1e-2)

    snr_ref = [-54.7270,  -53.5149,  -52.3028,  -51.0907,  -49.8786,  -48.6664,
               -47.4543,  -46.2422,  -45.0301,  -43.8180,  -42.6058,  -41.3937,
               -40.1816,  -38.9695,  -37.7574,  -36.5452,  -35.3331,  -34.1210,
               -32.9089,  -31.6967,  -30.4846,  -29.2725,  -28.0604,  -26.8483,
               -25.6361,  -24.4240,  -23.2119,  -21.9998,  -20.7877,  -19.5755,
               -18.3634,  -17.1513,  -15.9392,  -14.7270,  -13.5149,  -12.3028,
               -11.0907,   -9.8786,   -8.6664,   -7.4543,   -6.2422,   -5.0301,
               -3.8180,    -2.6058,   -1.3937,   -0.1816,    1.0305,    2.2426,
               +3.4548,     4.6669,    5.8790,    7.0911,    8.3032,    9.5154,
               +10.7275,   11.9396,   13.1517,   14.3638,   15.5759,   16.7881,
               +18.0002,   19.2123,   20.4244,   21.6365,   22.8485,   24.0606,
               +25.2727,   26.4847,   27.6967,   28.9087,   30.1206,   31.3324,
               +32.5442,   33.7558,   34.9673,   36.1785,   37.3895,   38.6002,
               +39.8103,   41.0198,   42.2285,   43.4360,   44.6421,   45.8462,
               +47.0478,   48.2458,   49.4393,   50.6266,   51.8058,   52.9741,
               +54.1277,   55.2617,   56.3694,   57.4405,   58.4607,   59.4074,
               +60.2442,   60.9031,   61.2360,   60.8103]
    amp_ref = [-120.0000, -118.7879, -117.5758, -116.3636, -115.1515, -113.9394,
               -112.7273, -111.5152, -110.3030, -109.0909, -107.8788, -106.6667,
               -105.4545, -104.2424, -103.0303, -101.8182, -100.6061,  -99.3939,
                -98.1818,  -96.9697,  -95.7576,  -94.5455,  -93.3333,  -92.1212,
                -90.9091,  -89.6970,  -88.4848,  -87.2727,  -86.0606,  -84.8485,
                -83.6364,  -82.4242,  -81.2121,  -80.0000,  -78.7879,  -77.5758,
                -76.3636,  -75.1515,  -73.9394,  -72.7273,  -71.5152,  -70.3030,
                -69.0909,  -67.8788,  -66.6667,  -65.4545,  -64.2424,  -63.0303,
                -61.8182,  -60.6061,  -59.3939,  -58.1818,  -56.9697,  -55.7576,
                -54.5455,  -53.3333,  -52.1212,  -50.9091,  -49.6970,  -48.4848,
                -47.2727,  -46.0606,  -44.8485,  -43.6364,  -42.4242,  -41.2121,
                -40.0000,  -38.7879,  -37.5758,  -36.3636,  -35.1515,  -33.9394,
                -32.7273,  -31.5152,  -30.3030,  -29.0909,  -27.8788,  -26.6667,
                -25.4545,  -24.2424,  -23.0303,  -21.8182,  -20.6061,  -19.3939,
                -18.1818,  -16.9697,  -15.7576,  -14.5455,  -13.3333,  -12.1212,
                -10.9091,   -9.6970,   -8.4848,   -7.2727,   -6.0606,   -4.8485,
                 -3.6364,   -2.4242,   -1.2121,         0]
    k0_ref = [3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6300,    3.6300,    3.6300,    3.6300,    3.6299,    3.6299,
              3.6298,    3.6298,    3.6297,    3.6295,    3.6293,    3.6291,
              3.6287,    3.6283,    3.6277,    3.6270,    3.6260,    3.6246,
              3.6228,    3.6205,    3.6173,    3.6131,    3.6075,    3.6000,
              3.5899,    3.5762,    3.5576,    3.5320,    3.4961,    3.4447,
              3.3690,    3.2518,    3.0562,    2.6817]
    k1_ref = [3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6301,    3.6301,    3.6301,
              3.6301,    3.6301,    3.6301,    3.6300,    3.6300,    3.6300,
              3.6300,    3.6299,    3.6299,    3.6298,    3.6297,    3.6296,
              3.6295,    3.6293,    3.6290,    3.6286,    3.6282,    3.6275,
              3.6267,    3.6256,    3.6242,    3.6223,    3.6198,    3.6164,
              3.6120,    3.6062,    3.5984,    3.5881,    3.5744,    3.5561,
              3.5318,    3.4993,    3.4557,    3.3969,    3.3171,    3.2075,
              3.0547,    2.8367,    2.5142,    2.0040]
    sigma_e2_ref = [0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,    0.3634,
                    0.3634,    0.3633,    0.3633,    0.3633,    0.3633,    0.3632,
                    0.3630,    0.3628,    0.3624,    0.3616,    0.3603,    0.3579,
                    0.3536,    0.3459,    0.3319,    0.3059]

    amp = np.linspace(-120, 0, 100)
    order = 4
    osr = 64
    f0 = 0.333
    Hinf = 1.2

    ntf = synthesizeNTF(order, osr, 2, Hinf, f0)
    snr_pred, amp_pred, k0, k1, sigma_e2 = predictSNR(ntf, osr, amp, f0)

    # Delete garbage data
    for check in (np.isinf, np.isnan):
        i = check(snr_ref)
        snr_ref = np.delete(snr_ref, np.where(i))
        amp_ref = np.delete(amp_ref, np.where(i))
        i = check(snr_pred)
        snr_pred = np.delete(snr_pred, np.where(i))
        amp_pred = np.delete(amp_pred, np.where(i))

    assert np.allclose(snr_pred, snr_ref, atol=1e-2, rtol=1e-3)
    assert np.allclose(amp_pred, amp_ref, atol=1e-2, rtol=5e-4)
    assert np.allclose(k0, k0_ref, atol=1e-3, rtol=1e-2)
    assert np.allclose(k1, k1_ref, atol=1e-3, rtol=50e-2)
    assert np.allclose(sigma_e2, sigma_e2_ref, atol=1e-3, rtol=1e-2)
