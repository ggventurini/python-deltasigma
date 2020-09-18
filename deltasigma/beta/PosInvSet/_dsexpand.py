# -*- coding: utf-8 -*-
# _dsexpand.py
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
from typing import Tuple, Optional

def dsexpand(s, c, k, n=None, o=None)->Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    function [sp, op] = dsexpand(s,c,k,n,o). Expand points s outward from c
    by a factor 1+k and update the offsets associated with normals n.
    """

    shift = np.tile([c[:, 0]], (1, np.shape(s)[1]))
    sp = shift + (1+k) * (s-shift)
    if (n is not None):
        op = o + k* (o+c*n)
    else:
        op = None

    return sp, op