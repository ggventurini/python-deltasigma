# -*- coding: utf-8 -*-
# _ds_f1f2.py
# This module provides the ds_f1f2 function.
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

"""This module provides the ds_f1f2() function.
"""

def ds_f1f2(OSR=64, f0=0., complex_flag=False):
    """[f1, f2] = ds_f1f2(OSR=64, f0=0, complex_flag=0)
    This function has no original docstring.
    """
    if complex_flag:
        f1 = f0 - 0.5/OSR
        f2 = f0 + 0.5/OSR
    else:
        if f0 > 0.25/OSR:
            f1 = f0 - 0.25/OSR
            f2 = f0 + 0.25/OSR
        else:
            f1 = 0.
            f2 = 0.5/OSR
    return f1, f2

def test_ds_f1f2():
    """Test function for ds_f1f2()"""
    f0 = 1e3
    OSR = 128
    cf = False
    t1f1, t1f2 = ds_f1f2(OSR, f0, cf)
    t2f1, t2f2 = ds_f1f2(OSR, f0, not cf)
    t3f1, t3f2 = ds_f1f2(OSR, .24/OSR, cf)
    assert (t1f1, t1f2) == (f0 - 0.25/OSR, f0 + 0.25/OSR)
    assert (t2f1, t2f2) == (f0 - 0.5/OSR, f0 + 0.5/OSR)
    assert (t3f1, t3f2) == (0., 0.5/OSR)

