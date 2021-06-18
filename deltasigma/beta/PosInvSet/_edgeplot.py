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
        plt.plot([s[0, :], s[0, 0]], [s[1, :], s[1, 0]], fmt)
        plt.grid()

    elif np.shape(s)[0] == 3:
        fig = plt.figure()
        tbl = [[1, 1, 2], 
               [4, 1, 3],
               [3, 2, 3],
               [2, 9, 9]]

        for p in range(3):
            ax = fig.add_subplot(2, 2, tbl[p, 0])
            x = tbl[p, 1]
            y = tbl[p, 2]
            ax.set_xlabel('x' + str(x))
            ax.set_ylabel('x' + str(y))

            for i in range(np.shape(e)[1]):
                p1 = s[:, e[0, i]]
                p2 = s[:, e[1, i]]
                ax.plot([p1[x], p2[x]], [p1[y], p2[y]], fmt)
            
            ax.grid(True)

        ax1 = fig.add_subplot(2, 2, tbl[3, 0], projection='3d')
        ax1.set_xlabel('x1')
        ax1.set_ylabel('x2')
        ax1.set_zlabel('x3')

        for i in range(np.shape(e)[1]):
            p1 = s[:, e[0, i]]
            p2 = s[:, e[1, i]]
            ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], fmt)

        ax1.grid(True)


    elif np.shape(s)[0] == 4:
        fig = plt.figure()
        tbl = [[1, 1, 2, 3],
               [2, 1, 2, 4],
               [3, 1, 3, 4],
               [4, 2, 3, 4]]

        for p in range(4):
            ax = fig.add_subplot(2, 2, tbl[p, 0], projection='3d')
            x = tbl[p, 1]
            y = tbl[p, 3]
            z = tbl[p, 4]
            ax.set_xlabel('x' + str(x))
            ax.set_ylabel('x' + str(y))
            ax.set_zlabel('x' + str(z))

            for i in range(np.shape(e)[1]):
                p1 = s[:, e[0, i]]
                p2 = s[:, e[1, i]]
                ax.plot([p1[x], p2[x]], [p1[y], p2[y]], [p1[z], p2[z]], fmt)

            ax.grid(True)
    else:
        pass
