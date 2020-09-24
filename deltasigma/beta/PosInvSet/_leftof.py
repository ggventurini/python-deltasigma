# -*- coding: utf-8 -*-
# _leftof.py
# Module providing the bilogplot function
# Copyright 2020 Yuki Fukuda
# This file is part of python-deltasigma(forked).
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.

import numpy as np

def leftof(p:np.ndarray, a:np.ndarray, b:np.ndarray)->bool:
    """function y=leftof(p,a,b)
    Return 1 if the point p is to the left of the line ab.
    For a n x 2 list of points p, return a vector of results.
    """

    p = p - np.tile(a[0, :], (1, np.shape(p)[0]))
    b = b - a
    y = b[0]*p[:, 1] > b[1]*p[:, 0]
    
    return y
    