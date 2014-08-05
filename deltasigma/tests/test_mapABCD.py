# -*- coding: utf-8 -*-
# test_mapABCD.py
# This module provides the tests for the mapABCD function.
# Copyright 2014 Giuseppe Venturini
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

"""This module provides the test class for the mapABCD() function.
"""

from __future__ import division

import unittest
import numpy as np

from nose.tools import raises

from deltasigma import mapABCD, realizeNTF, synthesizeNTF, stuffABCD

class TestMapABCD(unittest.TestCase):
    """Test class for mapABCD()"""

    def setUp(self):
        self.orders = (2, 3, 4, 5)
        self.osr = 32
        self.nlev = 2
        self.f0s = (0., 0.25)
        self.Hinf = 1.5
        self.forms = ('CRFB', 'CRFF', 'CIFB', 'CIFF', 'CRFFD', 'CRFBD', 'Stratos')
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
                   'CRFFD':{2:{'a':(0.7749, 0.2164),
                               'g':(0.),
                               'b':(1., 0., 1.),
                               'c':(1., 1.)
                              },
                            3:{'a':(0.7967, 0.2399, 0.0398),
                               'g':(0.0058),
                               'b':(1., 0., 0., 1.),
                               'c':(1., 1., 1.)
                              },
                            4:{'a':(0.8020, 0.2449, 0.0537, 0.0046),
                               'g':(0., 0.0069),
                               'b':(1., 0., 0., 0., 1.),
                               'c':(1., 1., 1., 1.)
                              },
                            5:{'a':(0.8023, 0.2443, 0.0570, 0.0071, 0.0002),
                               'g':(0.0028, 0.0079),
                               'b':(1., 0., 0., 0., 0., 1.),
                               'c':(1., 1., 1., 1., 1.)
                              }
                           },
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
                   'CRFFD':{2:{'a':(-1.2149e-15, -0.6667),
                               'g':(2.),
                               'b':(1., 0., 1.),
                               'c':(1., 1.)
                              },
                            4:{'a':(0.0000, -0.5585, -0.2164, -0.2164),
                               'g':(2., 2.),
                               'b':(1., 0., 0., 0., 1.),
                               'c':(1., 1., 1., 1.)
                              }
                           },
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
        self.res = res

    def test_mapABCD_1(self):
        """Test function for mapABCD() 1/2"""
        for f0 in self.f0s:
            for form in self.forms:
                for order in self.orders:
                    if f0 != 0. and order % 2 == 1:
                        # odd-order pass band modulator
                        continue
                    # Optimized zero placement
                    print("Testing form: %s, order: %d, f0: %f" % \
                          (form, order, f0))
                    ntf = synthesizeNTF(order, self.osr, 2, self.Hinf, f0)
                    a1, g1, b1, c1 = realizeNTF(ntf, form)
                    # we check realize NTF too
                    self.assertTrue(np.allclose(a1, self.res[f0][form][order]['a'],
                                    atol=1e-4, rtol=1e-3))
                    self.assertTrue(np.allclose(g1, self.res[f0][form][order]['g'],
                                    atol=1e-4, rtol=1e-3))
                    self.assertTrue(np.allclose(b1, self.res[f0][form][order]['b'],
                                    atol=1e-4, rtol=1e-3))
                    self.assertTrue(np.allclose(c1, self.res[f0][form][order]['c'],
                                    atol=1e-4, rtol=1e-3))
                    ABCD = stuffABCD(a1, g1, b1, c1, form)
                    a, g, b, c = mapABCD(ABCD, form)
                    self.assertTrue(np.allclose(a1, a, atol=1e-4, rtol=1e-3))
                    self.assertTrue(np.allclose(g1, g, atol=1e-4, rtol=1e-3))
                    self.assertTrue(np.allclose(b1, b, atol=1e-4, rtol=1e-3))
                    self.assertTrue(np.allclose(c1, c, atol=1e-4, rtol=1e-3))

    @raises(ValueError)
    def test_mapABCD_2(self):
        """Test function for mapABCD() 2/2"""
        ABCD = np.array([[0., 0., 0.91561444, -0.91561444],
                         [1., 0., 0., -1.42857142],
                         [0., 1., 0., 0.]])
        mapABCD(ABCD, 'DUMMY') 
