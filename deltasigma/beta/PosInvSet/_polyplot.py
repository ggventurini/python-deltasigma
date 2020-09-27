# -*- coding: utf-8 -*-
# _polyplot.py
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
import matplotlib.pyplot as plt

def polyplot(p, fmt = '-'):
    """
    function polyplot(p,fmt)
    Plot a polygon given by the point list p, with format fmt
    """

    n = np.shape(p)[1]

    if (np.array(p[:, n-1]) == np.array(p[:, 0])).all():
        plt.plot(p[0, :], p[1, :], fmt)
    else:
        plt.plot([p[0, :], p[0, 0]], [p[1, :], p[1, 0]], fmt)
