# -*- coding: utf-8 -*-
# _edgeplot.py
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
from mpl_toolkits.mplot3d import Axes3D

def edgeplot(e, s, fmt:str = '-'):
    if np.shape(s)[0] == 2:
        plt.plot([s[0, :] s[0, 0]], [s[1, :], s[1, 0]], fmt)
        plt.grid()

    elif np.shape(s)[0] == 3:
        fig = plt.figure()
        tbl = [[1, 1, 2], 
               [4, 1, 3],
               [3, 2, 3],
               [2, 9, 9]]

    elif np.shape(s)[0] == 4:

    else:
        pass