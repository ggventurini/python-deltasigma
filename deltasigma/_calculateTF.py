# -*- coding: utf-8 -*-
# _calculateTF.py
# Module providing the calculateTF function
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

"""Module providing the calculateTF() function
"""

import numpy as np
from scipy.signal import lti, ss2zpk
from ._constants import eps
from ._partitionABCD import partitionABCD
from ._utils import minreal

def calculateTF(ABCD, k=1.):
	"""Calculate the NTF and STF of a delta-sigma modulator given the loop filter
	ABCD matrix, assuming a quantizer gain of ``k``.

	**Parameters:**

	ABCD : array_like,
               The ABCD matrix that describes the system.

	k : float, optional
            The quantizer gain. If not specified, a default value of 1 is used.

	**Returns:**

	(NTF, STF) : a tuple of two LTI objects.

	**Example:**

	Realize a fifth-order modulator with the cascade-of-resonators structure,
	feedback form. Calculate the ABCD matrix of the loop filter and verify
	that the NTF and STF are correct.

	.. code-block:: python

	    from deltasigma import *
	    H = synthesizeNTF(5, 32, 1)
	    a, g, b, c = realizeNTF(H)
	    ABCD = stuffABCD(a,g,b,c)
	    ntf, stf = calculateTF(ABCD)

	From which we get:

	``H``::

	         (z -1) (z^2 -1.997z +1) (z^2 -1.992z +0.9999)      
	    --------------------------------------------------------
	     (z -0.7778) (z^2 -1.796z +0.8549) (z^2 -1.613z +0.665) 

	coefficients::

	    a: 0.0007, 0.0084, 0.055, 0.2443, 0.5579
	    g: 0.0028, 0.0079
	    b: 0.0007, 0.0084, 0.055, 0.2443, 0.5579, 1.0
	    c: 1.0, 1.0, 1.0, 1.0, 1.0

	ABCD matrix::

	    [[  1.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
	        0.00000000e+00   6.75559806e-04  -6.75559806e-04]
	     [  1.00000000e+00   1.00000000e+00  -2.79396240e-03   0.00000000e+00
	        0.00000000e+00   8.37752565e-03  -8.37752565e-03]
	     [  1.00000000e+00   1.00000000e+00   9.97206038e-01   0.00000000e+00
	        0.00000000e+00   6.33294166e-02  -6.33294166e-02]
	     [  0.00000000e+00   0.00000000e+00   1.00000000e+00   1.00000000e+00
	       -7.90937431e-03   2.44344030e-01  -2.44344030e-01]
	     [  0.00000000e+00   0.00000000e+00   1.00000000e+00   1.00000000e+00
	        9.92090626e-01   8.02273699e-01  -8.02273699e-01]
	     [  0.00000000e+00   0.00000000e+00   0.00000000e+00   0.00000000e+00
	        1.00000000e+00   1.00000000e+00   0.00000000e+00]]

	NTF::

	         (z -1) (z^2 -1.997z +1) (z^2 -1.992z +0.9999)      
	    --------------------------------------------------------
	     (z -0.7778) (z^2 -1.796z +0.8549) (z^2 -1.613z +0.665) 

	STF::

	    1

	"""
	A, B, C, D = partitionABCD(ABCD)
	if B.shape[1] > 1:
		B1 = B[:, 0]
		B2 = B[:, 1]
		B1 = B1.reshape((B1.shape[0], 1)) if len(B1.shape) == 1 else B1
		B2 = B2.reshape((B2.shape[0], 1)) if len(B2.shape) == 1 else B2
	else:
		B1 = B
		B2 = B

	# Find the noise transfer function by forming the closed-loop
	# system (sys_cl) in state-space form.
	Acl = A + k * np.dot(B2, C)
	Bcl = np.hstack((B1 + k*B2*D[0, 0], B2))
	Ccl = k*C
	Dcl = np.array((k*D[0, 0], 1.))
	Dcl = Dcl.reshape((1, Dcl.shape[0])) if len(Dcl.shape) == 1 else Dcl
	tol = min(1e-3, max(1e-6, eps**(1/ABCD.shape[0])))
	# input #0 is the signal
	# input #1 is the quantization noise
	stf_p, stf_z, stf_k  = ss2zpk(Acl, Bcl, Ccl, Dcl, input=0)
	ntf_p, ntf_z, ntf_k = ss2zpk(Acl, Bcl, Ccl, Dcl, input=1)
	stf = lti(stf_p, stf_z, stf_k)
	ntf = lti(ntf_p, ntf_z, ntf_k)
	stf_min, ntf_min = minreal((stf, ntf), tol)
	return ntf_min, stf_min

def test_calculateTF():
	"""Test function for calculateTF()"""
	from ._utils import cplxpair
	ABCD = [[1.000000000000000, 0., 0., 0.044408783846879, -0.044408783846879],
	        [0.999036450096481, 0.997109907515262, -0.005777399147297, 0., 0.499759089304780],
	        [0.499759089304780, 0.999036450096481, 0.997109907515262,  0., -0.260002096136488],
	        [0,                 0,                 1.000000000000000,  0, -0.796730400347216]]
	ABCD = np.array(ABCD)
	ntf, stf = calculateTF(ABCD)
	ntf_zeros, ntf_poles = np.roots(ntf.num), np.roots(ntf.den)
	stf_zeros, stf_poles = np.roots(stf.num), np.roots(stf.den)
	mntf_poles = np.array((1.498975311463384, 1.102565142679772, 0.132677264750882))
	mntf_zeros = np.array((0.997109907515262 + 0.075972576202904j, 
	                       0.997109907515262 - 0.075972576202904j, 
	                       1.000000000000000 + 0.000000000000000j)
	                      )
	mstf_zeros = np.array((-0.999999999999996,))
	mstf_poles = np.array((1.498975311463384, 1.102565142679772, 0.132677264750882))
	# for some reason, sometimes the zeros are in different order.
	ntf_zeros, mntf_zeros = cplxpair(ntf_zeros), cplxpair(mntf_zeros)
	stf_zeros, mstf_zeros = cplxpair(stf_zeros), cplxpair(mstf_zeros)
	ntf_poles, mntf_poles = cplxpair(ntf_poles), cplxpair(mntf_poles)
	stf_poles, mstf_poles = cplxpair(stf_poles), cplxpair(mstf_poles)
	assert np.allclose(ntf_zeros, mntf_zeros, rtol=1e-5, atol=1e-8)
	assert np.allclose(ntf_poles, mntf_poles, rtol=1e-5, atol=1e-8)
	assert np.allclose(stf_zeros, mstf_zeros, rtol=1e-5, atol=1e-8)
	assert np.allclose(stf_poles, mstf_poles, rtol=1e-5, atol=1e-8)
