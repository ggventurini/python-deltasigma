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
import timeit
from deltasigma import synthesizeNTF, realizeNTF, stuffABCD

# setup variables, wrap in functino for timeit

def setup_env():
    OSR = 32
    H = synthesizeNTF(5, OSR, 1)
    N = 8192
    f = 85
    u = 0.5*np.sin(2*np.pi*f/N*np.arange(N))
    a, g, b, c = realizeNTF(H, 'CRFB')
    ABCD = stuffABCD(a, g, b, c, form='CRFB')
    return u, H, ABCD

def test1():
    import _simulateDSM_numba
    u, H, ABCD = setup_env()
    _simulateDSM_numba.simulateDSM(u, H, 2, 0)

def test2():
    import _simulateDSM_numba
    u, H, ABCD = setup_env()
    _simulateDSM_numba.simulateDSM(u, ABCD, 2, 0)

def test3():
    import _simulateDSM_python
    u, H, ABCD = setup_env()
    _simulateDSM_python.simulateDSM(u, H, 2, 0)

def test4():
    import _simulateDSM_python
    u, H, ABCD = setup_env()
    _simulateDSM_python.simulateDSM(u, ABCD, 2, 0)

if __name__ == '__main__':
    print(timeit.timeit("test1()", setup="from __main__ import test1",
          number=1))
    print(timeit.timeit("test2()", setup="from __main__ import test2",
          number=1))
    print(timeit.timeit("test3()", setup="from __main__ import test3",
          number=1))
    print(timeit.timeit("test4()", setup="from __main__ import test4",
      number=1))
