# -*- coding: utf-8 -*-
# _realizeNTF.py
# Module providing realizeNTF
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

"""Module providing realizeNTF()
"""

from __future__ import division, print_function
from warnings import warn

import numpy as np
from scipy.signal import tf2zpk

from ._calculateTF import calculateTF
from ._evalRPoly import evalRPoly
from ._evalTF import evalTF
from ._stuffABCD import stuffABCD
from ._utils import carray, cplxpair, _get_zpk

def realizeNTF(ntf, form='CRFB', stf=None):
	"""Convert an NTF into coefficients for the desired structure.

	**Parameters:**

	ntf : a zpk transfer function
	    The modulator noise transfer function (NTF)

	form : string
	    A structure identifier.

	Supported structures are:
	
	* *"CRFB"*: Cascade of resonators, feedback form.

 	* *"CRFF"*: Cascade of resonators,  feedforward form.

	* *"CIFB"*: Cascade of integrators,  feedback form. 

	* *"CIFF"*: Cascade of integrators,  feedforward form.

 	* *"CRFBD"*: CRFB with delaying quantizer.

 	* *"PFF"*: Parallel feed-forward.

	* *"Stratos"*: A CIFF-like structure with non-delaying resonator feedbacks, contributed to the MATLAB delta sigma toolbox in 2007 by Jeff Gealow

	See the accompanying documentation (:ref:`topologies-diagrams`) for
	block diagrams of each structure.

	.. note:
	The order of the NTF zeros must be (real, complex conj. pairs).
	The order of the zeros is used when mapping the NTF onto the chosen topology.

	stf : a zpk transfer function, optional
	    the Signal Transfer Function
	
	**Returns:**

	a, g, b, c : tuple of ndarrays
	    the coefficients for the desired structure

	**Example:**

	Determine the coefficients for a 5th-order modulator with the
	cascade-of-resonators structure, feedback (CRFB) form.::

	    from deltasigma import synthesizeNTF, realizeNTF
	    H = synthesizeNTF(5, 32, 1)
	    a, g, b, c = realizeNTF(H,'CRFB')

	Returns the values::

	    a: 0.0007, 0.0084, 0.055, 0.2443, 0.5579
	    g: 0.0028, 0.0079
	    b: 0.0007, 0.0084, 0.055, 0.2443, 0.5579, 1.0
	    c: 1.0, 1.0, 1.0, 1.0, 1.0

	"""

	# The basic idea is to equate the loop filter at a set of
	# points in the z-plane to L1  =  1-1/ntf at those points.

	# Code common to all functions\
	if (hasattr(ntf, 'inputs') and not ntf.inputs == 1) or \
	   (hasattr(ntf, 'outputs') and not ntf.outputs == 1):
		raise TypeError("Only SISO transfer functions can be evaluated.")

	ntf_z, ntf_p, _k = _get_zpk(ntf)
	ntf_z = carray(ntf_z)
	ntf_p = carray(ntf_p)
	order = max(ntf_p.shape)
	order2 = int(np.floor(order/2))
	odd = (order - 2*order2 == 1)
	
	a = np.zeros((1, order))
	g = np.zeros((1, order2))
	b = np.zeros((1, order + 1))
	c = np.ones((1, order))
	T = np.zeros((order, order*2), dtype='complex')

	# instead of asuming enforce roots order
	ntf_z = cplxpair(ntf_z)[::-1]

	N = 200
	min_distance = 0.09
	zSet = np.zeros((N,), dtype='complex')
	j = 0
	for i in range(1, N + 1):
		z = 1.1 * np.exp(2j*np.pi*i/N)
		if np.all(np.abs(ntf_z - z) > min_distance):
			zSet[j] = z
			j = j + 1
	zSet = zSet[:j]
	if form == 'CRFB':
		# Find g
		# Assume the roots are ordered, real first, then cx conj. pairs 
		for i in range(0, order2):
			g[0, i] = 2*(1 - np.real(ntf_z[2*i + 1*odd]))
		L1 = np.zeros((1, order*2), dtype='complex')
		# Form the linear matrix equation a*T* = L1
		for i in range(0, 2*order):
			z = zSet[i]
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			Dfactor = (z - 1)/z
			product = 1
			for j in range(order, 0 + 1*odd, -2):
				product = z/evalRPoly(ntf_z[j - 2:j], z)*product
				T[j - 1, i] = product * Dfactor
				T[j - 2, i] = product
			if odd:
				T[0, i] = product/(z - 1)
		a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
		if stf is None:
			b[0, :order] = a[0, :]
			b[0, order] = 1
	elif form == 'CRFF':
		# Find g
		# Assume the roots are ordered, real first, then cx conj. pairs
		for i in range(0, order2):
			g[0, i] = 2*(1 - np.real(ntf_z[2*i + 1*odd]))
		L1 = np.zeros((1, order*2), dtype='complex')
		# Form the linear matrix equation a*T*=L1
		for i in range(0, 2*order):
			z = zSet[i]
			# L1(z) = 1-1/H(z)
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			if odd:
				Dfactor = z - 1
				product = 1/Dfactor
				T[0, i] = product
			else:
				Dfactor = (z - 1)/z
				product = 1.
			for j in range(0 + 1*odd, order, 2):
				product = z/evalRPoly(ntf_z[j:j + 2], z)*product
				T[j, i] = product * Dfactor
				T[j + 1, i] = product
		a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
		if stf is None:
			b = np.hstack((np.atleast_2d(1), 
			               np.zeros((1, order - 1)), 
			               np.atleast_2d(1)
			             ))
	elif form == 'CIFB':
		if any(abs(np.real(ntf_z) - 1) > 0.001):
			warn("The ntf's zeros have had their real parts set to one.")
		ntf_z = 1 + 1j*np.imag(ntf_z)
		for i in range(order2):
			g[0, i] = np.imag(ntf_z[2*i + 1*odd])**2
		L1 = np.zeros((1, order*2), dtype='complex')
		for i in range(order*2):
			z = zSet[i]
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			Dfactor = z - 1.
			product = 1.
			for j in range(order - 1, 1*odd, -2):
				product = product/evalRPoly(ntf_z[j-1:j+1], z)
				T[j, i] = product * Dfactor
				T[j - 1, i] = product
			if odd:
				T[0, i] = product/(z - 1.)
		a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
		if stf is None:
			b[0, :order] = a
			b[0, order] = 1.
	elif form == 'CIFF':
		if (abs(np.real(ntf_z) - 1) > 0.001).any():
			warn("The ntf's zeros have had their real parts set to one.")
		ntf_z = 1. + 1j*np.imag(ntf_z)
		for i in range(order2):
			g[0, i] = np.imag(ntf_z[2*i + 1*odd])**2
		L1 = np.zeros((1, order*2), dtype='complex')
		for i in range(order*2):
			z = zSet[i]
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			Dfactor = (z - 1)
			if odd:
				product = 1./(z - 1.)
				T[0, i] = product
			else:
				product = 1
			for j in range(0 + 1*odd, order - 1, 2):
				product = product/evalRPoly(ntf_z[j:j + 2], z)
				T[j, i] = product * Dfactor
				T[j + 1, i] = product
		a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
		if stf is None:
			b = np.zeros((1, order+1))
			b[0, 0] = 1.
			b[0, -1] = 1.
	elif form == 'CRFBD':
		for i in range(order2):
			g[0, i] = 2*(1. - np.real(ntf_z[2*i + 1*odd]))
		L1 = np.zeros((1, order*2), dtype='complex')
		for i in range(order*2):
			z = zSet[i]
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			Dfactor = (z - 1.)
			product = 1./z
			for j in range(order - 1, 0 + 1*odd, -2):
				product = z/evalRPoly(ntf_z[j - 1:j + 1], z) * product
				T[j, i] = product * Dfactor
				T[j - 1, i] = product
			if (odd):
				T[0, i] = product*z/(z - 1.)
		a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
		if stf is None:
			b = np.hstack((a, np.zeros((1, 1))))
			b[0, order] = 1.
	elif form == 'CRFFD':
		raise ValueError('I am sorry to inform you, the CRFFD is buggy in' +
		                 ' MATLAB and its Python port is equally buggy.\n' +
		                 'To this day, it is unclear what the original' +
		                 ' author meant with this code.\n'
		                 'If you know how to fix this, feel free to mail' +
		                 'me at ggventurini+GITHUB@gmail.com')
		# Find g
    		# Assume the roots are ordered, real first, then cx conj. pairs
		for i in range(order2):
			g[0, i] = 2*(1. -np.real(ntf_z[2*i + 1*odd]))
		# zL1 = z*(1-1/H(z))
		zL1 = np.atleast_2d(zSet * (1. - 1.0/evalTF(ntf, zSet)))
		# Form the linear matrix equation a*T*=zL1
		for i in range(order*2):
			z = zSet[i]
			if odd:
				Dfactor = (z - 1.)/z
				product = 1./Dfactor
				T[0, i] = product
			else:
				Dfactor = z - 1
				product = 1
			for j in range(0 + 1*odd, order, 2):
				product = z/evalRPoly(ntf_z[j:j + 2], z)*product
				T[j, i] = product * Dfactor
				T[j + 1, i] = product
		a = -np.real(np.linalg.lstsq(T.T, zL1.T)[0]).T
		if stf is None:
			b = np.hstack((np.atleast_2d(1), 
			               np.zeros(1, order - 1), 
			               np.atleast_2d(1))
			             )
	elif form == 'PFF':
		warn("Untested code accessed: realizeNTF('PFF')")
		# Find g
		# Assume the roots are ordered, real first, then cx conj. pairs
		# with the secondary zeros after the primary zeros
		for i in range(order2):
			g[0, i] = 2*(1. - np.real(ntf_z[2*i + 1*odd]))
		# Find the dividing line between the zeros
		theta0 = np.abs(np.angle(ntf_z[0]))
		# !! 0.5 radians is an arbitrary separation !!
		i = np.flatnonzero(np.abs(np.abs(np.angle(ntf_z)) - theta0) > 0.5)
		if i.shape[0]:
			order_1 = i[0] - 1
		else:
			order_1 = 0
		order_2 = order - order_1
		if i.shape[0] !=  order_2:
			raise ValueError('For the PFF form, the NTF zeros must ' +
			                 'be sorted into primary and secondary zeros')
		odd_1 = order_1 % 2
		odd_2 = order_2 % 2
		L1 = np.zeros((1, order*2), dtype='complex')
		# Form the linear matrix equation a*T*=L1
		for i in range(order*2):
			z = zSet[i]
			# L1(z) = 1-1/H(z)
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			if odd_1:
				Dfactor = z - 1
				product = 1./Dfactor
				T[0, i] = product
			else:
				Dfactor = (z - 1)/z
				product = 1
			for j in range(odd_1, order_1, 2):
				product = z/evalRPoly(ntf_z[j:j + 2], z) * product
				T[j, i] = product * Dfactor
				T[j + 1, i] = product
			if odd_2:
				Dfactor = z - 1
				product = 1./Dfactor
				T[order_1, i] = product
			else:
				Dfactor = (z - 1)/z
				product = 1.
			for j in range(order_1 + odd_2, order, 2):
				product = z/evalRPoly(ntf_z[j:j + 2], z)*product
				T[j, i] = product * Dfactor
				T[j + 1, i] = product
		a = -np.real(np.linalg.lstsq(T.T, zL1.T)[0]).T
		if stf is None:
			b = np.hstack((np.atleast_2d(1), 
			               np.zeros((1, order_1 - 1)), 
			               np.atleast_2d(1), 
			               np.zeros(1, order_2 - 1), 
			               np.atleast_2d(1)
			             ))
	elif form == 'Stratos':
		# code copied from case 'CRFF':
		# Find g
		# Assume the roots are ordered, real first, then cx conj. pairs
		for i in range(order2):
			g[0, i] = 2*(1 - np.real(ntf_z[2*i + 1*odd]))
		# code copied from case 'CIFF':
		L1 = np.zeros((1, order*2), dtype='complex')
		# Form the linear matrix equation a*T*=L1
		for i in range(order*2):
			z = zSet[i]
			# L1(z) = 1-1/H(z)
			L1[0, i] = 1. - evalRPoly(ntf_p, z)/evalRPoly(ntf_z, z)
			Dfactor = (z - 1.)
			if odd:
				product = 1./(z - 1.)
				T[0, i] = product
			else:
				product = 1.
			for j in range(0 + 1*odd, order - 1, 2):
				product = product/evalRPoly(ntf_z[j:j + 2], z)
				T[j, i] = product * Dfactor
				T[j + 1, i] = product
		a = -np.real(np.linalg.lstsq(T.T, L1.T)[0]).T
		if stf is None:
			b = np.hstack((
			               np.atleast_2d(1),
			               np.zeros((1, order - 1)), 
			               np.atleast_2d(1)
			             ))
	if not stf is None:
		# Compute the TF from each feed-in to the output 
		# and solve for coefficients which yield the best match
		# THIS CODE IS NOT OPTIMAL,  in terms of computational efficiency.
		stfList = []
		for i in range(0, order + 1):
			bi = np.zeros(1, order + 1)
			bi[i] = 1
			ABCD = stuffABCD(a, g, bi, c, form)
			if form[2:4] == 'FF':
				ABCD[0, order + 1] = - 1
			_, stfListi = calculateTF(ABCD)
			stfList.append(stfListi)
		# Build the matrix equation b A  =  x and solve it.
		A = np.zeros((order + 1, max(zSet.shape)))
		for i in range(0, order + 1):
			A[i, :] = evalTF(stfList[i], zSet)
		x = evalTF(stf, zSet)
		x = x[:].T
		b = np.real(x/A)
	return a.squeeze(), g.squeeze(), b.squeeze(), c.squeeze()
	
def test_realizeNTF():
	"""Test function for realizeNTF()"""
	from ._synthesizeNTF import synthesizeNTF
	orders = (2, 3, 4, 5)
	osr = 32
	nlev = 2
	f0s = (0., 0.25)
	Hinf = 1.5
	forms = ('CRFB', 'CRFF', 'CIFB', 'CIFF', 'CRFBD', 'Stratos')
	res = {0.:{'CIFB':{2:{'a':(0.2164, 0.7749),
	                      'g':(0, ),
	                      'b':(0.2164, 0.7749, 1.0000),
	                      'c':(1., 1. )
	                     },
	                   3:{'a':(0.0444, 0.2843, 0.8025),
	                      'g':(0.0058, ),
	                      'b':(0.0444, 0.2843, 0.8025, 1.),
	                      'c':(1., 1., 1.)
	                     },
	                   4:{'a':(0.0062, 0.0655, 0.3042, 0.8089),
	                      'g':(0., 0.0069),
	                      'b':(0.0062, 0.0655, 0.3042, 0.8089, 1.),
	                      'c':(1., 1., 1., 1)
	                     },
	                   5:{'a':(0.0007, 0.0095, 0.0731, 0.3100, 0.81309),
	                      'g':(0.0028, 0.0079),
	                      'b':(0.0007, 0.0095, 0.0731, 0.3100, 0.8130, 1.),
	                      'c':(1., 1., 1., 1., 1.)
	                     }
	                  },
	           'CRFB':{2:{'a':(0.2164, 0.5585),
	                      'g':(0, ),
	                      'b':(0.2164, 0.5585, 1.0000),
	                      'c':(1., 1. )
	                     },
	                   3:{'a':(0.0444, 0.2399, 0.5569),
	                      'g':(0.0058, ),
	                      'b':(0.0444, 0.2399, 0.5569, 1.),
	                      'c':(1., 1., 1.)
	                     },
	                   4:{'a':(0.0062, 0.0530, 0.2449, 0.5571),
	                      'g':(0, 0.0069),
	                      'b':(0.0062, 0.0530, 0.2449, 0.5571, 1.),
	                      'c':(1., 1., 1., 1)
	                     },
	                   5:{'a':(0.0007, 0.0084, 0.0550, 0.2443, 0.5579),
	                      'g':(0.0028, 0.0079),
	                      'b':(0.0007, 0.0084, 0.0550, 0.2443, 0.5579, 1.),
	                      'c':(1., 1., 1., 1., 1.)
	                     }
	                   },
	           'CRFF':{2:{'a':(0.5585, 0.2164),
	                      'g':(0, ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. )
	                     },
	                   3:{'a':(0.5569, 0.2399, 0.0412),
	                      'g':(0.0058, ),
	                      'b':(1., 0., 0., 1.),
	                      'c':(1., 1., 1.)
	                     },
	                   4:{'a':(0.5571, 0.2449, 0.0492, 0.0046),
	                      'g':(0, 0.0069),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1)
	                     },
	                   5:{'a':(0.5579, 0.2443, 0.0505, 0.0071, 0.0003),
	                      'g':(0.0028, 0.0079),
	                      'b':(1., 0., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1., 1.)
	                     }
	                   },
	           'CIFF':{2:{'a':(0.7749, 0.2164),
	                      'g':(0., ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. )
	                     },
	                   3:{'a':(0.8025, 0.2843, 0.0398),
	                      'g':(0.0058, ),
	                      'b':(1., 0., 0., 1.),
	                      'c':(1., 1., 1.)
	                     },
	                   4:{'a':(0.8089, 0.3042, 0.0599, 0.0041),
	                      'g':(0., 0.0069),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1)
	                     },
	                   5:{'a':(0.8130, 0.3100, 0.0667, 0.0080, 0.0001),
	                      'g':(0.0028, 0.0079),
	                      'b':(1., 0., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1., 1.)
	                     }
	                  },
	           'CRFBD':{2:{'a':(0.2164, 0.7749),
	                       'g':(0, ),
	                       'b':(0.2164, 0.7749, 1.),
	                       'c':(1., 1. )
	                      },
	                    3:{'a':(0.0444, 0.2399, 0.7967),
	                       'g':(0.0058, ),
	                       'b':(0.0444, 0.2399, 0.7967, 1.),
	                       'c':(1., 1., 1.)
	                      },
	                    4:{'a':(0.0062, 0.0592, 0.2449, 0.8020),
	                       'g':(0, 0.0069),
	                       'b':(0.0062, 0.0592, 0.2449, 0.8020, 1.),
	                       'c':(1., 1., 1., 1)
	                      },
	                    5:{'a':(0.0007, 0.0084, 0.0633, 0.2443, 0.8023),
	                       'g':(0.0028, 0.0079),
	                       'b':(0.0007, 0.0084, 0.0633, 0.2443, 0.8023, 1.),
	                       'c':(1., 1., 1., 1., 1.)
	                      }
	                   },
	           #'CRFFD':{2:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         3:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         4:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         5:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           }
	           #        }
	           'Stratos':{2:{'a':(0.7749, 0.2164),
	                       'g':(0, ),
	                       'b':(1, 0., 1.),
	                       'c':(1., 1. )
	                      },
	                    3:{'a':(0.7967, 0.2796, 0.0398),
	                       'g':(0.0058, ),
	                       'b':(1., 0., 0., 1.),
	                       'c':(1., 1., 1.)
	                      },
	                    4:{'a':(0.8020, 0.2987, 0.0579, 0.0042),
	                       'g':(0, 0.0069),
	                       'b':(1., 0., 0., 0., 1.),
	                       'c':(1., 1., 1., 1)
	                      },
	                    5:{'a':(0.8023, 0.3013, 0.0643, 0.0075, 0.0001),
	                       'g':(0.0028, 0.0079),
	                       'b':(1., 0., 0., 0., 0., 1.),
	                       'c':(1., 1., 1., 1., 1.)
	                      }
	                   }
	          },
	      0.25:{'CIFB':{2:{'a':(0.3333, 2.0000),
	                      'g':(1., ),
	                      'b':(0.3333, 2., 1.0000),
	                      'c':(1., 1. )
	                     },
	                   4:{'a':(-3.5585, 2.4503, 5.22512, 4.0000),
	                      'g':(1., 1.),
	                      'b':(-3.5585, 2.4503, 5.22512, 4.0000, 1.),
	                      'c':(1., 1., 1., 1)
	                     }
	                  },
	           'CRFB':{2:{'a':(-0.6667, 0.6667),
	                      'g':(2.0, ),
	                      'b':(-0.6667, 0.6667, 1.0000),
	                      'c':(1., 1. )
	                     },
	                   4:{'a':(-0.2164, 0., -0.5585, 0.5585),
	                      'g':(2.0, 2.0),
	                      'b':(-0.2164, 0.0, -0.5585, 0.5585, 1.),
	                      'c':(1., 1., 1., 1)
	                     }
	                   },
	           'CRFF':{2:{'a':(0.6667, -0.6667),
	                      'g':(2.0, ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. )
	                     },
	                   4:{'a':(0.5585, -0.5585, 0., -0.2164),
	                      'g':(2., 2.),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1)
	                     }
	                   },
	           'CIFF':{2:{'a':(2., 0.3333),
	                      'g':(1., ),
	                      'b':(1., 0., 1.),
	                      'c':(1., 1. )
	                     },
	                   4:{'a':(4., 5.2251, 2.4503, -3.5585),
	                      'g':(1., 1.),
	                      'b':(1., 0., 0., 0., 1.),
	                      'c':(1., 1., 1., 1)
	                     }
	                  },
	           'CRFBD':{2:{'a':(-0.6667, 0.),
	                       'g':(2.0, ),
	                       'b':(-0.6667, 0., 1.),
	                       'c':(1., 1. )
	                      },
	                    4:{'a':(-0.2164, -0.2164, -0.5585, 0.),
	                       'g':(2., 2.),
	                       'b':(-0.2164, -0.2164, -0.5585, 0., 1.),
	                       'c':(1., 1., 1., 1)
	                      }
	                   },
	           #'CRFFD':{2:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           },
	           #         4:{'a':(),
	           #            'g':(),
	           #            'b':(),
	           #            'c':()
	           #           }
	           #        }
	           'Stratos':{2:{'a':(0., -0.6667),
	                       'g':(2., ),
	                       'b':(1, 0., 1.),
	                       'c':(1., 1. )
	                      },
	                    4:{'a':(0., -0.7749, 0., 0.2164),
	                       'g':(2.0, 2.0),
	                       'b':(1., 0., 0., 0., 1.),
	                       'c':(1., 1., 1., 1)
	                      }
	                   }
	          }
	      }

	for f0 in f0s:
		for form in forms:
			for order in orders:
				if f0 != 0. and order % 2 == 1:
					# odd-order pass band modulator
					continue
				# Optimized zero placement
				print("Testing form: %s, order: %d, f0: %f" % \
				      (form, order, f0))
				ntf = synthesizeNTF(order, osr, 2, Hinf, f0)
				a, g, b, c = realizeNTF(ntf, form)
				print(a, g, b, c)
				assert np.allclose(a, res[f0][form][order]['a'], 
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(g, res[f0][form][order]['g'], 
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(b, res[f0][form][order]['b'], 
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(c, res[f0][form][order]['c'], 
				            atol=1e-4, rtol=1e-3)
	return 
