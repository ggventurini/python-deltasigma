# -*- coding: utf-8 -*-
# _outsideConvex.py
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

def outsideConvex(x:np.ndarray, n:np.ndarray, o:np.ndarray, tol:float=0)->np.ndarray:
    if (np.shape(x)[1] * np.shape(n)[1] < 100000):
        i = np.sign(np.sum(np.array(n).T*x + np.tile(o[0, :], (np.shape(x)[1], 1)).T > True))
    else:
        chunk = np.round(100000/np.shape(x)[1])
        i = np.zeros(1, np.shape(x)[1])

        for j in range(0, np.shape(x)[1], chunk):
            j2 = np.min(j+chunk-1, np.shape(x)[1])
            i[j:j2] = np.sign(np.sum(np.array(n).T * x[:, j:j2] + np.tile(o[0, :], (j2-j+1, 1)).T > tol))