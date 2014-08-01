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

def predictSNR(ntf, OSR=64, amp=None, f0=0.):
    """Predict the SNR curve of a binary delta-sigma modulator.

    The prediction is performed using the describing function method of Ardalan
    and Paulos [2]_ .

    **Parameters:**
    
    ntf : lti object, or zpk or (num, den) or (A,B,C,D) tuples
        The noise transfer function specifying the modulator.
    
    OSR : scalar, optional
        The oversampling ratio, defaults to 64.
    
    amp : ndarray-like, optional
        The magnitudes to be used for the input signal. They are expressed in 
        dB, where 0 dB means a full-scale (peak value = 1) sine wave. 
        Defaults to [-120 -110 ... -20 -15 -10 -9 -8 ... 0].
    
    f0 : scalar, optional 
        The normalized input signal frequency. Defaults to 0.
    
    **Notes:**
    
    The band of interest is defined by the oversampling ratio (``OSR``) and
    the center frequency (``f0``).

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
        band_of_interest = np.linspace(0, np.pi/OSR, Nb)
    else:
        band_of_interest = np.linspace(2*np.pi*(f0 - 0.25/OSR), 2*np.pi*(f0 + 0.25/OSR), Nb)
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
        snr[0, n] = dbp(0.5*u[n]**2/(np.sum(h**2)/(OSR*Nb)*sigma_e2[0, n]))
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

