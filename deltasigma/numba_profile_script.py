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
import cProfile
import pstats

def test_script():
    OSR = 32
    H = synthesizeNTF(5, OSR, 1)
    N = 8192
    f = 85
    u = 0.5*np.sin(2*np.pi*f/N*np.arange(N))
    a, g, b, c = realizeNTF(H, 'CRFB')
    ABCD = stuffABCD(a, g, b, c, form='CRFB')

    _simulateDSM_numba.simulateDSM(u, H, 2, 0)

        #_simulateDSM_numba.simulateDSM(u, ABCD, 2, 0)

        #_simulateDSM_python.simulateDSM(u, H, 2, 0)

        #_simulateDSM_python.simulateDSM(u, ABCD, 2, 0)

cProfile.run('test_script()', 'numba_stats')
p = pstats.Stats('numba_stats')
p.sort_stats('time').print_stats(25)
