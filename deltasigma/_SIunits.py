# -*- coding: utf-8 -*-
# _SIunits.py
# This module provides the SIunits function.
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

"""This module provides the SIunits() function, which, given an input,
returns factor (think "engineering notation") and an alphabetic
suffix.
"""

from __future__ import division
import numpy as np 
from ._constants import eps

def SIunits(x):
    """Calculates the factor for representing ``x`` in engineering 
    notation.

    The factors and suffixes supported are:

    +----------+---+---------+--------+---+---------+
    |   1e-3   | m |  milli  |   1e3  | k |   kilo  |
    +----------+---+---------+--------+---+---------+
    |   1e-6   | u |  micro  |   1e6  | M |   mega  |
    +----------+---+---------+--------+---+---------+
    |   1e-9   | n |   nano  |   1e9  | G |   giga  |
    +----------+---+---------+--------+---+---------+
    |   1e-12  | p |   pico  |  1e12  | T |   tera  |
    +----------+---+---------+--------+---+---------+
    |   1e-15  | f |  femto  |  1e15  | P |   peta  |
    +----------+---+---------+--------+---+---------+
    |   1e-18  | a |   atto  |  1e18  | E |    exa  |
    +----------+---+---------+--------+---+---------+
    |   1e-21  | z |  zepto  |  1e21  | Z |   zeta  |
    +----------+---+---------+--------+---+---------+
    |   1e-24  | y |  yocto  |  1e24  | Y |  yotta  |
    +----------+---+---------+--------+---+---------+

    **Parameters:**

    x : scalar or array_like
            The number for which the engineering notation factor and suffix
        are to be calculated.

    **Returns:**

    factor : float or list of floats
             the engineering notation factor(s)

    prefix : string or list of strings
             the engineering notation unit prefix(es)

    **Example**::

        a = 3300.
        unit = 'g'
        f, p = SIunits(a)
        print "Float 'a' in engineering notation: %.3f %s%s" % (a/f, p, unit)

    Prints::

        Float 'a' in engineering notation: 3.300 kg
    """
    prefixes_n = ('', 'm', 'u', 'n', 'p', 'f', 'a', 'z', 'y')
    prefixes_p = ('', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    scalar_input = np.isscalar(x)
    if scalar_input:
        x = [x]
    factor = np.ones(len(x))
    prefix = ['']*len(x)
    for i in range(len(x)):
        if x[i] != 0:
            p = int(np.floor(np.log10(abs(x[i]))/3 + eps))
            if p >= 0:
                p = min(p, len(prefixes_p))
                prefix[i] = prefixes_p[p]
                factor[i] = 10**(3*p)
            elif p < 0:
                p = min(-p, len(prefixes_n)-1)
                prefix[i] = prefixes_n[p]
                factor[i] = 10**(-3*p)
        else:
            factor[i], prefix[i] = 0, ''
    if scalar_input:
        factor, prefix = factor[0], prefix[0]
    return factor, prefix

def test():
    """Test function for SIunits()"""
    tv = (0, 1, 1e3, 2100312.24, .32545, 21e-9, 34e-12, 9569300e-12)
    correct = (0, ''), (1, ''), (1e3, 'k'), (1e6, 'M'), (1e-3, 'm'), (1e-9, 'n'), (1e-12, 'p'), (1e-6, 'u')
    f, p = SIunits(tv)
    res = zip(f,p)
    for r, c in zip(res, correct):
        assert r[0] == c[0] and r[1] == c[1]
    # test scalars
    tv = 2100312.24
    correct = (1e6, 'M')
    f, p = SIunits(tv)
    assert f == correct[0] and p == correct[1]

