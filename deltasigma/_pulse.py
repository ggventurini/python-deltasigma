# -*- coding: utf-8 -*-
# _pulse.py
# This module provides the pulse function.
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

"""This module provides the pulse() function, which calculates the sampled
pulse response of a CT system.
"""

from __future__ import division
import collections
import numpy as np
from scipy.signal import step2, lti

from ._utils import lcm, rat
from ._utils import _is_zpk, _is_num_den, _is_A_B_C_D

def pulse(S, tp=(0., 1.), dt=1., tfinal=10., nosum=False):
    """Calculate the sampled pulse response of a CT system.

    ``tp`` may be an array of pulse timings, one for each input, or even a
    simple 2-elements tuple.

    **Parameters:**

    S : sequence
        A sequence of LTI objects specifying the system.

    The sequence S should be assembled so that ``S[i][j]`` returns the
    LTI system description from input ``i`` to the output ``j``.

    In the case of a MISO system, a unidimensional sequence ``S[i]``
    is also acceptable.

    tp : array-like
        An (n, 2) array of pulse timings

    dt : scalar
        The time increment

    tfinal : scalar
        The time of the last desired sample

    nosum : bool
        A flag indicating that the responses are not to be summed

    **Returns:**

    y : ndarray
        The pulse response
    """
    tp = np.asarray(tp)
    if len(tp.shape) == 1:
        if not tp.shape[0] == 2:
            raise ValueError("tp is not (n, 2)-shaped")
        tp = tp.reshape((1, tp.shape[0]))
    if len(tp.shape) == 2:
        if not tp.shape[1] == 2:
            raise ValueError("tp is not (n, 2)-shaped")

    # Compute the time increment
    dd = 1
    for tpi in np.nditer(tp.T.copy(order='C')):
        _, di = rat(tpi, 1e-3)
        dd = lcm(di, dd)

    _, ddt = rat(dt, 1e-3)
    _, df = rat(tfinal, 1e-3)
    delta_t = 1./lcm(dd, lcm(ddt, df))
    delta_t = max(1e-3, delta_t)    # Put a lower limit on delta_t
    if (isinstance(S, collections.Iterable) and len(S)) \
       and (isinstance(S[0], collections.Iterable) and len(S[0])) \
           and (isinstance(S[0][0], lti) or _is_zpk(S[0][0]) or _is_num_den(S[0][0]) \
                or _is_A_B_C_D(S[0][0])):
        pass
    else:
        S = list(zip(S)) #S[input][output]
    y1 = None
    for Si in S:
        y2 = None
        for So in Si:
            _, y2i = step2(So, T=np.arange(0., tfinal + delta_t, delta_t))
            if y2 is None:
                y2 = y2i.reshape((y2i.shape[0], 1, 1))
            else:
                y2 = np.concatenate((y2,
                                     y2i.reshape((y2i.shape[0], 1, 1))),
                                    axis=1)
        if y1 is None:
            y1 = y2
        else:
            y1 = np.concatenate((y1, y2), axis=2)

    nd = int(np.round(dt/delta_t, 0))
    nf = int(np.round(tfinal/delta_t, 0))
    ndac = tp.shape[0]

    ni = len(S) # number of inputs

    if ni % ndac != 0:
        raise ValueError('The number of inputs must be divisible by the number of dac timings.')
        # Original comment from the MATLAB sources:
        # This requirement comes from the complex case, where the number of inputs
        # is 2 times the number of dac timings. I think this could be tidied up.

    # nis: Number of inputs grouped together with a common DAC timing
    # (2 for the complex case)
    nis = int(ni/ndac)

    # notice len(S[0]) is the number of outputs for us
    if not nosum: # Sum the responses due to each input set
        y = np.zeros((np.ceil(tfinal/float(dt)) + 1, len(S[0]), nis))
    else:
        y = np.zeros((np.ceil(tfinal/float(dt)) + 1, len(S[0]), ni))

    for i in range(ndac):
        n1 = int(np.round(tp[i, 0]/delta_t, 0))
        n2 = int(np.round(tp[i, 1]/delta_t, 0))
        z1 = (n1, y1.shape[1], nis)
        z2 = (n2, y1.shape[1], nis)
        yy = + np.concatenate((np.zeros(z1), y1[:nf-n1+1, :, i*nis:(i + 1)*nis]), axis=0) \
             - np.concatenate((np.zeros(z2), y1[:nf-n2+1, :, i*nis:(i + 1)*nis]), axis=0)
        yy = yy[::nd, :, :]
        if not nosum: # Sum the responses due to each input set
            y = y + yy
        else:
            y[:, :, i] = yy.reshape(yy.shape[0:2])
    return y

