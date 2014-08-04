# -*- coding: utf-8 -*-
# test_utils.py
# This module provides the tests for the assorted functions in the utils
# module.
# Copyright 2014 Giuseppe Venturini
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

"""This module provides the test class for the assorted functions in the
_utils module.
"""

import unittest
import numpy as np
import deltasigma as ds

from deltasigma import synthesizeNTF

from deltasigma._utils import *
from deltasigma._utils import _get_zpk, _get_num_den, _cell_like_list
from deltasigma._constants import eps

def test_rat():
    """Test function for rat()"""
    import numpy.random as rnd
    for _ in range(10):
        n, d = rnd.randint(1, 5000), rnd.randint(1, 5000)
        fr = float(n) / float(d)
        nt, dt = rat(fr, tol=200e-6)
        assert np.allclose(
            (n / float(d) - nt / float(dt),), (0.,), atol=200e-6, rtol=1e-12)


def test_gcd_lcm():
    """Test function for gcd() and lcm"""
    a, b = 36, 31721
    tlcm, tgcd = 1141956, 1
    assert lcm(a, b) == tlcm
    assert gcd(a, b) == tgcd


def test_mfloor():
    """Test function for mfloor()"""
    tv = np.linspace(-1, 1, 10)
    tres = np.zeros(tv.shape)
    tres[tv >= 0] = np.floor(tv[tv >= 0])
    tres[tv < 0] = -np.ceil(np.abs(tv[tv < 0]))
    tresf = mfloor(tv)
    assert np.allclose(tres, tresf, atol=1e-8, rtol=1e-5)
    # w complex values.
    tv = [-1.9,  -0.2,  3.4,  5.6,  7.0,  2.4 + 3.6j]
    tres = [-2.0, -1.0, 3.0, 5.0, 7.0, (2 + 3j)]
    assert mfloor(tv) == tres


def test_carray():
    """Test function for carray()"""
    # don't touch 1d arrays
    a = np.arange(10)
    b = a.copy()
    assert np.all(carray(a) == b)
    # reshape 0-d arrays
    a = np.array(1)
    b = a.reshape((-1,))
    assert np.all(carray(a) == b)
    # convert lists and tuples to 1d arrays
    a = np.arange(10)
    c = a.tolist()
    b = tuple(c)
    assert np.all(a == carray(b))
    assert np.all(a == carray(c))
    # convert scalars to 1d arrays
    a = np.array((1,))
    b = 1
    assert np.all(carray(b) == a)


def test_cplxpair():
    """Test function for cplxpair()
    """
    a = np.array(
        [1 + eps * 20j, 1.1 + 2j, 1.1 - (2 + 50 * eps) * 1j, .1 + (1 + 99 * eps) * .2j, .1 - .2j])
    assert np.allclose(
        cplxpair(a), np.array(
            [0.1 - 0.2j, 0.1 + 0.2j, 1.1 - 2.j, 1.1 + 2.j, 1.0 + 0.j]),
        atol=100 * eps)
    # Mismatched imaginary part with correct sign
    a = (1, 1+2j, 1-2.4j) # this should fail. The complex sing. are not conj.
    try:
        cplxpair(a)
        # ValueError should be raised!
        assert False
    except ValueError:
        assert True
    # Mismatched real part with correct imaginary part
    a = (1, 1+2j, 2-2j) # this should fail. The complex sing. are not conj.
    try:
        cplxpair(a)
        # ValueError should be raised!
        assert False
    except ValueError:
        assert True
    # Mismatched real part with mismatched sign imaginary part
    a = (1, 1+2j, 2+2j) # this should fail. The complex sing. are not conj.
    try:
        cplxpair(a)
        # ValueError should be raised!
        assert False
    except ValueError:
        assert True

def test_diagonal_indices():
    """Test function for diagonal_indices()
    """
    a = np.arange(1, 26)
    a = a.reshape((5, 5))
    d = a[diagonal_indices(a)]
    dp1 = a[diagonal_indices(a, 1)]
    dm1 = a[diagonal_indices(a, -1)]
    assert np.allclose(d, (1, 7, 13, 19, 25))
    assert np.allclose(dp1, (2, 8, 14, 20))
    assert np.allclose(dm1, (6, 12, 18, 24))
    a = np.arange(1, 21)
    a = a.reshape((4, 5))
    d = a[diagonal_indices(a)]
    dp1 = a[diagonal_indices(a, 1)]
    dm1 = a[diagonal_indices(a, -1)]
    assert np.allclose(d, (1, 7, 13, 19))
    assert np.allclose(dp1, (2, 8, 14, 20))
    assert np.allclose(dm1, (6, 12, 18))
    a = a.reshape((5, 4))
    d = a[diagonal_indices(a)]
    dp1 = a[diagonal_indices(a, 1)]
    dm1 = a[diagonal_indices(a, -1)]
    assert np.allclose(d, (1, 6, 11, 16))
    assert np.allclose(dp1, (2, 7, 12))
    assert np.allclose(dm1, (5, 10, 15, 20))


def test_circshift():
    """Test function for circshift()
    """
    A = np.arange(1, 10).reshape((3, 3))
    shift = np.array((1, -1))
    Ashifted = circshift(A, shift)
    Ares = np.array(((8, 9, 7), (2, 3, 1), (5, 6, 4)))
    assert np.allclose(Ashifted, Ares)
    a = np.array([1, 2, 3, 4, 5]) # test for shift being a scalar
    assert np.allclose(circshift(a, 1), a[np.array((4, 0, 1, 2, 3), )])


def test_save_input_form():
    """Test function for save_input_form() and restore_input_form()"""
    a1 = np.arange(1, 26)
    a2 = a1.reshape((1, -1))
    form_a1 = save_input_form(a1)
    assert a1.shape == restore_input_form(a2, form_a1).shape
    assert np.allclose(a1, restore_input_form(a2, form_a1))
    b1 = 5
    b2 = carray(b1)
    form_b1 = save_input_form(b1)
    assert np.isscalar(restore_input_form(b2, form_b1))
    assert b1 == restore_input_form(b2, form_b1)


def test_get_zpk():
    """Test function for _get_zpk()"""
    from scipy.signal import zpk2ss
    z = (.4,)
    p = (.9, .1 + .2j, .1 - .2j)
    k = .4
    zt, pt, kt = _get_zpk(zpk2ss(z, p, k))
    assert np.allclose(zt, z, atol=1e-8, rtol=1e-5)
    assert np.allclose(pt, p, atol=1e-8, rtol=1e-5)
    assert np.allclose((kt, ), (k, ), atol=1e-8, rtol=1e-5)
    zt, pt, kt = _get_zpk(zpk2tf(z, p, k))
    assert np.allclose(zt, z, atol=1e-8, rtol=1e-5)
    assert np.allclose(pt, p, atol=1e-8, rtol=1e-5)
    assert np.allclose((kt, ), (k, ), atol=1e-8, rtol=1e-5)
    zt, pt, kt = _get_zpk(lti(z, p, k))
    assert np.allclose(zt, z, atol=1e-8, rtol=1e-5)
    assert np.allclose(pt, p, atol=1e-8, rtol=1e-5)
    assert np.allclose((kt, ), (k, ), atol=1e-8, rtol=1e-5)


def test_get_num_den():
    """Test function for _get_num_den()"""
    from scipy.signal import zpk2ss
    z = (.4,)
    p = (.9, .1 + .2j, .1 - .2j)
    k = .4
    num, den = zpk2tf(z, p, k)
    numt, dent = _get_num_den((num, den))  # num, den
    assert np.allclose(numt, num, atol=1e-8, rtol=1e-5)
    assert np.allclose(dent, den, atol=1e-8, rtol=1e-5)
    numt, dent = _get_num_den((z, p, k))  # zpk
    assert np.allclose(numt, num, atol=1e-8, rtol=1e-5)
    assert np.allclose(dent, den, atol=1e-8, rtol=1e-5)
    numt, dent = _get_num_den(lti(z, p, k))  # LTI
    assert np.allclose(numt, num, atol=1e-8, rtol=1e-5)
    assert np.allclose(dent, den, atol=1e-8, rtol=1e-5)
    assert len(numt.shape) == 1
    assert len(dent.shape) == 1
    A, B, C, D = zpk2ss(z, p, k)  # A,B,C,D
    D = np.atleast_2d(D)
    numt, dent = _get_num_den((A, B, C, D))  # A,B,C,D
    assert np.allclose(numt, num, atol=1e-8, rtol=1e-5)
    assert np.allclose(dent, den, atol=1e-8, rtol=1e-5)
    ABCD = np.vstack((
                     np.hstack((A, B)),
                     np.hstack((C, D))
                     ))
    numt, dent = _get_num_den(ABCD)  # A,B,C,D
    assert np.allclose(numt, num, atol=1e-8, rtol=1e-5)
    assert np.allclose(dent, den, atol=1e-8, rtol=1e-5)
    H = [[1],[2]] # check no 0-length arrays are returned
    numt, dent = _get_num_den(H)
    assert np.allclose(numt, 1, atol=1e-8, rtol=1e-5)
    assert np.allclose(dent, 2, atol=1e-8, rtol=1e-5)
    assert len(numt.shape) == 1
    assert len(dent.shape) == 1


def test_minreal():
    """Test function for minreal()"""
    from scipy.signal import zpk2ss
    z = np.array([0, 1])
    p = np.array([1, .5 + .5j, .5 - .5j])
    k = 1
    zt = z[:-1]
    pt = p[1:]
    pt.sort()
    inobj = (z, p, k)
    ltiobj = minreal(inobj)
    zeros, poles, gain = _get_zpk(ltiobj)
    poles.sort()
    assert np.allclose(zeros, zt, atol=1e-8, rtol=1e-5)
    assert np.allclose(poles, poles, atol=1e-8, rtol=1e-5)
    assert np.allclose((k,), (gain,), atol=1e-8, rtol=1e-5)
    inobj = zpk2tf(z, p, k)
    ltiobj = minreal(inobj)
    zeros, poles, gain = _get_zpk(ltiobj)
    poles.sort()
    assert np.allclose(zeros, zt, atol=1e-8, rtol=1e-5)
    assert np.allclose(poles, poles, atol=1e-8, rtol=1e-5)
    assert np.allclose((k,), (gain,), atol=1e-8, rtol=1e-5)
    inobj = zpk2ss(z, p, k)
    ltiobj = minreal(inobj)
    zeros, poles, gain = _get_zpk(ltiobj)
    poles.sort()
    assert np.allclose(zeros, zt, atol=1e-8, rtol=1e-5)
    assert np.allclose(poles, poles, atol=1e-8, rtol=1e-5)
    assert np.allclose((k,), (gain,), atol=1e-8, rtol=1e-5)
    inobj = lti(z, p, k)
    ltiobj = minreal(inobj)
    zeros, poles, gain = _get_zpk(ltiobj)
    poles.sort()
    assert np.allclose(zeros, zt, atol=1e-8, rtol=1e-5)
    assert np.allclose(poles, poles, atol=1e-8, rtol=1e-5)
    assert np.allclose((k,), (gain,), atol=1e-8, rtol=1e-5)


def test_pretty_lti():
    """Test function for pretty_lti()"""
    H = synthesizeNTF(5, 32, 1)
    tv = pretty_lti(H)
    res = '      (z^2 - 1.992z + 0.9999) (z^2 - 1.997z + 1) (z - 1)     \n' + \
          '-------------------------------------------------------------\n' + \
          ' (z^2 - 1.613z + 0.665) (z^2 - 1.796z + 0.8549) (z - 0.7778) '
    assert res == tv
    assert int(pretty_lti(([], [], 2))) == 2
    assert pretty_lti([[1], [], 2]) == '2 (z - 1)'
    assert pretty_lti([[], [.22222222], 2]) == \
        '      2       \n--------------\n (z - 0.2222) '
    assert pretty_lti(((0, 0, 1), (1, 2-1j, 2+1j, 2-1j, 2+1j), 5)) == \
        '         z^2 (z - 1)        \n' + \
        '5 --------------------------\n' + \
        '   (z^2 - 4z + 5)^2 (z - 1) '


def test_cell_like_list():
    """Test function for _cell_like_list()"""
    res = [[[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]],
           [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]],
           [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]]]
    assert _cell_like_list((3, 4, 5), 1) == res


def test_mround():
    """Test function for mfloor()"""
    tv = np.linspace(-1, 1, 21)
    tres = np.zeros(tv.shape)
    tres[:6] = -1
    tres[-6:] = 1
    tresm = mround(tv)
    assert np.allclose(tres, tresm, atol=1e-8, rtol=1e-5)
    # w complex values.
    tv = [-1.9, -0.5, -0.2, 3.4, 4.5, 5.6, 7.0, 2.4 + 3.6j]
    tres = [-2.0, -1.0,  0.0, 3.0, 5.0, 6.0, 7.0, 2.0 + 4.0j]
    assert mround(tv) == tres
