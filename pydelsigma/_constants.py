# -*- coding: utf-8 -*-
# _constants.py
# This module holds a few constants that are built-in in MATLAB and that were
# not found in Python.
# Copyright 2013 Giuseppe Venturini
# This file is part of python-deltasigma.
#
# python-deltasigma is a 1:1 Python replacement of Richard Schreier's 
# MATLAB delta sigma toolbox (aka "delsigma"), upon which it is heavily based.
# The delta sigma toolbox is (c) 2009, Richard Schreier.
#
# python-deltasigma is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# LICENSE file for the licensing terms.

"""This module holds a few constants that are built-in in MATLAB and that were
not found in Python.
"""

import numpy as np

eps = np.finfo(np.float).eps # x86 2.22044604925e-16

