# -*- coding: utf-8 -*-
# _mapABCD.py
# Module providing the mapABCD function
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

"""Module providing the mapABCD() function
"""

from __future__ import division, print_function
import numpy as np

from ._utils import diagonal_indices

def mapABCD(ABCD, form='CRFB'):
    """Compute the coefficients for the specified structure.

    It is assumed that the ABCD matrix fits the topology.

    **Parameters:**

    ABCD : ndarray
        A state-space description of the modulator loop filter.
    form : str, optional
        See :func:`realizeNTF` for a list of supported structures.

    **Returns:**
    
    a : ndarray
        Feedback/feedforward coefficients from/to the quantizer. Length :math:`n`.
    g : ndarray
        Resonator coefficients. Length :math:`floor(n/2)`.
    b : ndarray
        Feed-in coefficients from the modulator input to each integrator. Length :math:`n + 1`.
    c : ndarray
        Integrator inter-stage coefficients. Length :math:`n`.

    .. seealso::

        * :func:`realizeNTF` for a list of supported structures.

        * :func:`stuffABCD`, the inverse function.

    """
    ABCD = np.copy(ABCD)
    order = ABCD.shape[0] - 1
    odd = order % 2
    even = 1 - odd
    diagonal = diagonal_indices(ABCD)
    subdiag = diagonal_indices(ABCD, -1)
    supdiag = [a[odd:order - 1:2] for a in diagonal_indices(ABCD, +1)]
    if form in ('CRFB', 'CIFB', 'CRFBD'):
        c = ABCD[subdiag]
        g = -ABCD[supdiag]
        if form == 'CRFB':
            dly = np.arange(1 + odd, order, 2)
            ABCD[dly, :] = ABCD[dly, :] \
                           - np.dot(np.diag(c[dly - 1]), ABCD[dly - 1, :])
        elif form == 'CRFBD':
            dly = np.arange(odd, order, 2)
            ABCD[dly, :] = ABCD[dly, :] + np.dot(np.diag(g), ABCD[dly + 1, :])
            if order > 2:
                coupl = np.arange(1 + even, order, 2)
                ABCD[coupl, :] = ABCD[coupl, :] \
                                 - np.dot(np.diag(c[coupl - 1]), ABCD[coupl - 1, :])
        a = -ABCD[:order, order + 1].T
        b = ABCD[:, order].T
    elif form == 'CRFF':
        a = np.zeros((order,))
        c = np.concatenate((
                            np.array((-ABCD[0, order + 1],)), 
                            ABCD[subdiag][:-1]
                          ))
        g = -ABCD[supdiag]
        if even:
            multg = np.arange(0, order, 2)
            ABCD[multg, :] = ABCD[multg, :] + np.dot(np.diag(g), ABCD[multg + 1, :])
        multc = np.arange(2, order, 2)
        ABCD[multc, :] = ABCD[multc, :] - np.dot(np.diag(c[multc]), ABCD[multc - 1, :])
        a[1:order:2] = ABCD[order, 1:order:2]
        for i in range(1, order, 2):
            ABCD[order, :] = ABCD[order, :] - a[i]*ABCD[i, :]
        a[:order:2] = ABCD[order, :order:2]
        b = ABCD[:, order].T
    elif form == 'CRFFD':
        #order=order - 1
        #odd=rem(order,2)
        #even=not  odd
        diagonal = diagonal_indices(ABCD[:order, :order])
        subdiag = diagonal_indices(ABCD[:order, :order], -1)
        supdiag = [a[1+odd:order:2] for a in 
                      diagonal_indices(ABCD[:order, :order], -1)]
        g = -ABCD[supdiag]
        c = np.vstack((-ABCD[0, order + 2], ABCD[subdiag]))
        a = np.zeros((1, order))
        for i in range(0, order, 2):
            a[i] = ABCD[order, i]
            ABCD[order, :] = ABCD[order, :] - a[i]*ABCD[i, :]
        a[1:order:2] = ABCD[order, 1:order:2]
        b = ABCD[:order + 1, order + 1].T
        for i in range(1, order, 2):
            b[i] = b[i] - c[i]*b[i - 1]
            if odd:
                b[i] = b[i] + g[(i - 1)/2]*b[i + 1]
        yscale = ABCD[order + 1, order]
        a = a*yscale
        b[-1] = b[-1]*yscale
    elif form in ('CIFF', 'Stratos'):
        a = ABCD[order, :order]
        c = np.concatenate((np.atleast_1d(-ABCD[0, order + 1]),
                            ABCD[subdiag][:-1]))
        g = -ABCD[supdiag]
        b = ABCD[:, order].T
    else:
        raise ValueError('Form %s is not yet supported.' % form)

    return a.squeeze(), g.squeeze(), b.squeeze(), c.squeeze()

def test_mapABCD():
	"""Test function for mapABCD()"""
	from ._realizeNTF import realizeNTF
	from ._synthesizeNTF import synthesizeNTF
	from ._stuffABCD import stuffABCD
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
				a1, g1, b1, c1 = realizeNTF(ntf, form)
				ABCD = stuffABCD(a1, g1, b1, c1, form)
				a, g, b, c = mapABCD(ABCD, form)
				assert np.allclose(a, res[f0][form][order]['a'],
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(g, res[f0][form][order]['g'],
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(b, res[f0][form][order]['b'],
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(c, res[f0][form][order]['c'],
				            atol=1e-4, rtol=1e-3)
				assert np.allclose(a1, a, atol=1e-4, rtol=1e-3)
				assert np.allclose(g1, g, atol=1e-4, rtol=1e-3)
				assert np.allclose(b1, b, atol=1e-4, rtol=1e-3)
				assert np.allclose(c1, c, atol=1e-4, rtol=1e-3)
	return
