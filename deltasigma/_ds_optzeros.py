# -*- coding: utf-8 -*-
# _ds_optzeros.py
# Module providing the optzeros function
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

"""Module providing the ds_optzeros() function
"""

import numpy as np
from numpy import sqrt

_oznopt = {
           1:{1:np.array((0.,)), 2:np.array((0.,))},
           2:{1:np.sqrt(np.array((1./3,))), 2:np.array((0.,))},
           3:{1:np.array((sqrt(3./5), 0.)), 2:np.array((sqrt(3./5), 0.))},
           4:{1:sqrt(np.array((3./7 + sqrt(9./49 - 3./35), 3./7 - sqrt(9./49 - 3./35)))),
              2:np.array((0., sqrt(5./7)))},
           5:{1:sqrt(np.array((5./9 + sqrt(25./81 - 5./21), 5./9 - sqrt(25./81 - 5./21), 0.))),
              2:sqrt(np.array((5./9 + sqrt(25./81 - 5./21), 5./9 - sqrt(25./81 - 5./21), 0.)))},
           6:{1:np.array((0.23862059, 0.66120988, 0.9324696)),
              2:sqrt(np.array((0., 7./11 + sqrt(56.)/33, 7./11 - sqrt(56.)/33)))},
           7:{1:np.array((0., 0.40584371, 0.74153078, 0.94910785)),
              2:np.array((0., 0.40584371, 0.74153078, 0.94910785))},
           8:{1:np.array((0.18343709, 0.52553345, 0.79666684, 0.96028993)),
              2:np.array((0., 0.50563161, 0.79017286, 0.95914731))},
           9:{1:np.array((0., 0.32425101, 0.61337056, 0.83603082, 0.9681602)),
              2:np.array((0., 0.32425101, 0.61337056, 0.83603082, 0.9681602))},
           10:{1:np.array((0.1834370913, 0.5255334458, 0.7966668433, 0.9602899327)),
              2:np.array((0., 0.41572267, 0.67208682, 0.86238894, 0.97342121))},
           11:{1:np.array((0., 0.26953955, 0.51909468, 0.73015137, 0.88706238, 0.97822864)),
               2:np.array((0., 0.26953955, 0.51909468, 0.73015137, 0.88706238, 0.97822864))},
           12:{1:np.array((0.12523875, 0.36783403, 0.58731921, 0.7699033, 0.90411753, 0.9815607)),
               2:np.array((0., 0.35222363, 0.58006251, 0.76647993, 0.90281326, 0.98132047))},
           13:{1:np.array((0., 0.23045331, 0.44849063, 0.64234828, 0.8015776, 0.91759824, 0.98418306)),
               2:np.array((0., 0.23045331, 0.44849063, 0.64234828, 0.8015776, 0.91759824, 
                           0.98418306))},
           14:{1:np.array((0.10806212, 0.31911586, 0.51525046, 0.68729392, 0.82720185, 
                       0.92843513, 0.98628389)),
               2:np.array((0., 0.30524384, 0.50836649, 0.6836066, 0.82537239, 
                       0.92772336, 0.98615167))}
           }

def ds_optzeros(n, opt=1):
    """A helper function for :func:`synthesizeNTF`

    Returns the zeros which minimize the in-band noise power of 
    a delta-sigma modulator's NTF.

    This function is not intended for direct use, but it is available
    for compliance with the Matlab Toolbox interface.

    **Parameters:**

    n : int
        The order of the modulator

    opt : int
        A flag which selects the kind of optimization to be employed
        for the zeros. A description of the possible values can be found
        in the doc for :func:`synthesizeNTF`.

    **Returns:**

    zeros : 1d-ndarray
        An array with the location of the zeros in the z plane, according
        to the specified optimization.

    """
    opt = int(opt)    
    if opt == 0:
        optZeros = np.zeros((np.ceil(n/2.), ))
    else:
        optZeros = _oznopt[n][opt]
    
    # Sort the zeros and replicate them.
    z = np.sort(optZeros)
    optZeros = np.zeros((n,))
    m = 0
    if n % 2 == 1:
        optZeros[0] = z[0]
        z = z[1:]
        m += 1
    for i in range(z.shape[0]):
        optZeros[m]     =  z[i]
        optZeros[m + 1] = -z[i]
        m += 2
    return optZeros

