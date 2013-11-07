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
from scipy.signal import zpk2tf, dimpulse, freqz
from scipy.special import erfinv
from scipy.interpolate import interp1d

from ._dbp import dbp

def predictSNR(ntf, R=64, amp=None, f0=0.):
    """
    Predict the SNR curve of a binary delta-sigma modulator by using the describing function 
    method of Ardalan and Paulos.
    
    Parameters:
    ===========
    
    ntf: lti or tf object, or zpk or (num, den) tuples
         The noise transfer function specifying the modulator.
    
    R: scalar, default 64, oversampling ratio
    
    amp: ndarray-like, default [-120 -110...-20 -15 -10 -9 -8 ... 0], the magnitudes to be used 
    for the input signal. They are expressed in dB, where 0 dB means a full-scale 
    (peak value = 1) sine wave. 
    
    f0, scalar, default 0 is the (notmalized) input signal frequency.
    
    Notes:
    ======
    
    The band of interest is defined by the oversampling ratio (R) and the center frequency (f0).

    The algorithm assumes that the amp vector is sorted in increasing order;
    once instability is detected, the remaining SNR values are set to -Inf.
    
    Output:
    =======

    snr: ndarray, a vector of SNR values (in dB)

    amp: ndarray, a vector of amplitudes (in dB)

    k0: ndarray, the quantizer signal gain

    k1: ndarray, the quantizer noise gain

    sigma_e2: scalar, the power of the quantizer noise (not in dB)
    

    Implementation details
    ======================

    The describing function method of A&P assumes that the quantizer processes
    signal and noise components separately. The quantizer is modelled as two
    (not necessarily equal) linear gains, k0 and k1, and an additive white
    gaussian noise source of power sigma_e2. k0, k1 and sigma_e2 are calculated
    as functions of the input.
    
    The modulator's loop filter is assumed to have nearly infinite gain at
    the test frequency.

    Future versions may accommodate STFs.
    """

    # extract num, den
    if (hasattr(ntf, 'inputs') and not ntf.inputs == 1) or \
       (hasattr(ntf, 'outputs') and not ntf.outputs == 1):
        raise TypeError("The supplied TF isn't a SISO transfer function.")
    if hasattr(ntf, 'num') and hasattr(ntf, 'den'):
        filt = hasattr(ntf, '__class__') and ntf.__class__.__name__ == 'TransferFunction'
        num = ntf.num[0][0] if filt else ntf.num
        den = ntf.den[0][0] if filt else ntf.den
    elif (hasattr(ntf, 'zeros') and hasattr(ntf, 'poles')) or \
         (hasattr(ntf, 'zero') and hasattr(ntf, 'pole')):
        # LTI objects have poles and zeros, 
        # TransferFunction-s have pole() and zero()
        zeros = ntf.zeros if hasattr(ntf, 'zeros') else ntf.zero()
        poles = ntf.poles if hasattr(ntf, 'poles') else ntf.pole()
        num, den = zpk2tf(zeros, poles, 1)
    elif hasattr(ntf, 'form') and ntf.form == 'coeff':
        num, den = ntf.num, ntf.den
    elif hasattr(ntf, 'form'):
        raise ValueError('%s: Unknown form: %s' % (__name__, ntf.form))
    elif hasattr(ntf, '__len__'):
        if len(ntf) == 3: # z, p, k
            zeros, poles, _ = ntf
            num, den = zpk2tf(zeros, poles, 1)
        elif len(ntf) == 2: # num, den
            num, den = ntf[0], ntf[1]
        else:
            raise TypeError('%s: Unknown transfer function %s' % (__name__, str(ntf)))
    else:
        raise TypeError('%s: Unknown transfer function %s' % (__name__, str(ntf)))

    Nb = 100
    if f0 == 0:
        band_of_interest = np.linspace(0, np.pi/R, Nb)
    else:
        band_of_interest = np.linspace(2*np.pi*(f0 - 0.25/R), 2*np.pi*(f0 + 0.25/R), Nb)
        XTAB = np.linspace(-2, 0, 21)
        YTAB = np.array([0.465759605169, 0.673669993877, 0.479046523571, 0.684266507626, 
                         0.493162959814, 0.695279479027, 0.508173644543, 0.70673173666, 
                         0.52414894104, 0.718647658825, 0.541165232658, 0.731052994728, 
                         0.559305548668, 0.743975520134, 0.578660130501, 0.757444560528, 
                         0.599327206612, 0.771491587162, 0.621413528919, 0.786150157452, 
                         0.645035266876, 0.801456093788, 0.670318901539, 0.817447543144, 
                         0.697402179241, 0.834165394306, 0.72643494606, 0.851653397083, 
                         0.757580637932, 0.869958162308, 0.791017174721, 0.889129817486, 
                         0.826938569546, 0.90922164917, 0.865556240082, 0.930291116238, 
                         0.907100915909, 0.952399373054, 0.951824009418, 0.975612223148, 
                         1.0, 1.0]).reshape(1,-1)

    if amp is None:
        amp = np.concatenate((
                              np.arange(- 120, -20 + 1, 10),
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
                m0 = interp1d(XTAB, YTAB[:, 1], kind='cubic')
                f = rho*m0(-rho**2) - K
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
            m1 = interp1d(XTAB, YTAB[:, 0], kind='cubic')
            sigma_e2[0, n] = 1 - u[n]**2/2 - 2/np.pi*m1(-rho**2)**2
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
            k1sigma1 = np.sqrt(2/np.pi) * m1(-rho**2)
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
