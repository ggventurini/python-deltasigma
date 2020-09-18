# -*- coding: utf-8 -*-
# _dotplot.py
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

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def dotplot(p:np.ndarray, fmt:str = '.'):
    
    if np.shape(p)[0] == 2:
        plt.plot(p[0, :], p[1, :], fmt)
    elif np.shape(p)[0] == 3:
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax1.plot(p[0, :], p[1, :], fmt)
        ax2 = fig.add_subplot(2, 2, 2, projection='3d')
        ax2.plot(p[0, :], p[1, :], p[2, :], fmt)
        ax3 = fig.add_subplot(2, 2, 4)
        ax3.plot(p[0, :], p[2, :], fmt)
        ax4 = fig.add_subplot(2, 2, 3)
        ax4.plot(p[1, :], p[2, :], fmt)
    elif np.shape(p)[0] == 4:
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1, projection='3d')
        ax1.plot(p[0, :], p[1, :], p[2, :], fmt)
        ax2 = fig.add_subplot(2, 2, 2, projection='3d')
        ax2.plot(p[0, :], p[1, :], p[3, :], fmt)
        ax3 = fig.add_subplot(2, 2, 3, projection='3d')
        ax3.plot(p[0, :], p[2, :], p[3, :], fmt)
        ax4 = fig.add_subplot(2, 2, 4, projection='3d')
        ax4.plot(p[1, :], p[2, :], p[3, :], fmt)
    else:
        pass