# -*- coding: utf-8 -*-
# _frespF1.py
# Module providing the frespF1 function
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

"""Module providing the frespF1() function
"""

from __future__ import division
import numpy as np
import pylab as plt

from ._dbv import dbv

def frespF1(f1, f=None, phi=1, plot=False):
    """Plot/calculate the frequency response of the F1 filter 
    in a Saramaki HBF at the points given in the optional f (n by 1) vector.
    """
    if f is None:
        f = np.linspace(0, 0.5)
    cos_w = np.cos(2*np.pi*f)
    F1 = 0.5
    for i in range(max(f1.shape)):
        F1 = F1 + f1[i] * ((cos_w/phi)**(2*i + 1))
    if plot:
        plt.plot(f, dbv(F1))
        plt.grid('on')
    fresp = F1
    return fresp
