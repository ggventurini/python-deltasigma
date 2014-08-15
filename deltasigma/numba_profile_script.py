# -*- coding: utf-8 -*-
# numba_timing_tests.py
# This script is to help time changes to _simulateDSM_numba.py
# Copyright 2014 Giuseppe Venturini & Shayne Hodge
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

"""This script is to help time changes to _simulateDSM_numba.py
"""

import numpy as np
from deltasigma import synthesizeNTF, realizeNTF, stuffABCD
import _simulateDSM_numba
import _simulateDSM_python
import cProfile
import pstats

def setup():
    OSR = 32
    H = synthesizeNTF(5, OSR, 1)
    N = 8192
    f = 85
    u = 0.5*np.sin(2*np.pi*f/N*np.arange(N))
    a, g, b, c = realizeNTF(H, 'CRFB')
    ABCD = stuffABCD(a, g, b, c, form='CRFB')
    return OSR, H, N, f, u, ABCD

def test_script_numba():
    OSR, H, N, f, u, ABCD = setup()
    #_simulateDSM_numba.simulateDSM(u, H, 2, 0)
    _simulateDSM_numba.simulateDSM(u, ABCD, 2, 0)

def test_script_cpython():
    OSR, H, N, f, u, ABCD = setup()
    #_simulateDSM_python.simulateDSM(u, H, 2, 0)
    _simulateDSM_python.simulateDSM(u, ABCD, 2, 0)

cProfile.run('test_script_numba()', 'numba_stats')
pn = pstats.Stats('numba_stats')
print('Numba stats')
pn.sort_stats('time').print_stats(10)
pn.sort_stats('cumtime').print_stats(10)

cProfile.run('test_script_cpython()', 'cpython_stats')
pc = pstats.Stats('cpython_stats')
print('CPython stats')
pc.sort_stats('time').print_stats(10)
