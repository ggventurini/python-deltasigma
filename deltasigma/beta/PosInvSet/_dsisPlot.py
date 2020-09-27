# -*- coding: utf-8 -*-
# _dsisPlot.py
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
from ._dotplot import dotplot
from ._polyplot import polyplot
from ._edgeplot import edgeplot

def dsisPlot(dbg, itn, order, x, s, e, ns, out):
    if dbg == 0:
        return
    elif dbg == 2:
        print("Iteration {}: {}/{} image vertices outside".format(itn, np.int(np.sum(out)), np.shape(ns)[1]))
        return
    else:
        pass

    if itn > 0:
        str_title = "Iteration {}: {}/{} image points outside".format(itn, np.int(np.sum(out)), np.shape(ns)[1])
    else:
        str_title = "Final object: {} vertices and {} edges".format(np.shape(s)[1], np.shape(e)[1])

    plt.close()

    if order == 2:
        dotplot(x)
        plt.grid()
        polyplot(s, 'k')
        dotplot(ns, ko)
        dotplot(ns[:, np.nonzero(out)[0]], 'rs')
        plt.title(str_title)
    elif order >= 3:
        dotplot(x)
        edgeplot(e, s)
        dotplot(ns, '+')
        dotplot(ns[:, np.nonzero(out)[0]], 'rs')
        plt.subplot(222)
        plt.title(str_title)
    else:
        pass
