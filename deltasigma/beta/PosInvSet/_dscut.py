# -*- coding: utf-8 -*-
# _dscut.py
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

def dscut(p1: np.ndarray, y1: np.ndarray, p2: np.ndarray, y2: np.ndarray)->np.ndarray:
    """Return the point p as which y = 0, assuming y varies linearly

    Parameters
    ----------


    Returns
    -------

    """

    k1 = y2 / (y2 - y1)
    k2 = 1 - k1
    n = np.shape(p1)[0]
    p = np.tile([k1[0, :]], (n, 1)) * p1 + np.tile([k2[0, :]], (n, 1)) * p2

    return p
