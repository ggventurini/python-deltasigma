# -*- coding: utf-8 -*-
# _uvar.py
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

def uvar(u, N)->np.ndarray:
    """
    un = uvar(u,N);	Compute a bounded sequence that has 90% of its 
    values at the extremes.
    """

    u1 = u[0]
    u2 = u[1]
    un = np.tile(u1, N).astype(np.float64)

    r = np.random.rand(N)

    ri2_tmp = (r > 0.55)
    ri2 = ri2_tmp * 1
    ri1 = ((r > 0.45) * 1) * ((~ri2_tmp) * 1)
    
    un[np.nonzero(ri1)[0]] = u1 + (u2 - u1) * 10 * (r[np.nonzero(ri1)[0]] - 0.45)
    un[np.nonzero(ri2)[0]] = u2

    return un
