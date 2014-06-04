# -*- coding: utf-8 -*-
# _utils.py
# Miscellaneous functions and stdlib wrappers for MATLAB functions
# that do not find a direct replacement in numpy/scipy.
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

""" Miscellaneous functions and wrappers for MATLAB functions
 that do not find a direct replacement in numpy/scipy.
"""

import collections
import fractions
from fractions import Fraction as Fr
import numpy as np
from scipy.signal import lti, ss2tf, ss2zpk, zpk2tf

from._constants import eps
from ._partitionABCD import partitionABCD


def rat(x, tol):
    """Rational fraction approximation.

    Calculate A and B such that:

    .. math::

      x = \\frac{A}{B} + \\epsilon

    where:

    .. math::

        |\\epsilon| < tol


    .. note:: A, B are of type 'int'
    """
    return Fr(float(x)).limit_denominator(int(1 / float(tol))).numerator, \
        Fr(float(x)).limit_denominator(int(1 / float(tol))).denominator

gcd = fractions.gcd

lcm = lambda a, b: int(a * b / float(gcd(a, b)))
lcm.__doc__ = """Calculate the Least Common Multiple of ``a`` and ``b``.\n"""


class empty(object):

    """An empty function used to hold attributes"""

    def __init__(self):
        pass


def mfloor(x):
    """Round ``x`` towards -Inf.

    This is a MATLAB-compatible floor function, numpy's ``floor()``
    behaves differently.

    If the elements of ``x`` are complex, real and imaginary parts are
    rounded separately.
    """
    iform = save_input_form(x)
    x = carray(x)

    def _mfloor(z):
        """Base function to generate the ufunc floor"""
        z = np.real_if_close(z)
        if np.iscomplex(z):
            return _mfloor(np.real(z)) + 1j * _mfloor(np.imag(z))
        return np.floor(z) if np.sign(z) >= 0 else -np.ceil(-z)
    _internal = np.frompyfunc(_mfloor, 1, 1)
    xf = np.array(_internal(x), dtype=x.dtype)
    return restore_input_form(xf, iform)


def carray(x):
    """Check that x is an ndarray. If not, try to convert it to ndarray.
    """
    if not hasattr(x, 'shape'):
        if not hasattr(x, '__len__'):
            x = np.array((x,))
        else:
            x = np.array(x)
    elif not len(x.shape):
        x = x.reshape((1,))
    else:
        pass  # nothing to do here
    return x


def cplxpair(x, tol=100):
    """
    Sort complex numbers into complex conjugate pairs.

    This function replaces MATLAB's cplxpair for vectors.
    """
    x = carray(x)
    x = np.atleast_1d(x.squeeze())
    x = x.tolist()
    x = [np.real_if_close(i, tol) for i in x]
    xreal = np.array(list(filter(np.isreal, x)))
    xcomplex = np.array(list(filter(np.iscomplex, x)))
    xreal = np.sort_complex(xreal)
    xcomplex = np.sort_complex(xcomplex)
    xcomplex_ipos = xcomplex[xcomplex.imag > 0.]
    xcomplex_ineg = xcomplex[xcomplex.imag <= 0.]
    if len(xcomplex_ipos) != len(xcomplex_ineg):
        raise ValueError("Complex numbers can't be paired.")
    res = []
    for i, j in zip(xcomplex_ipos, xcomplex_ineg):
        if not abs(i - np.conj(j)) < tol * eps:
            raise ValueError("Complex numbers can't be paired.")
        res += [j, i]
    return np.hstack((np.array(res), xreal))


def minreal(tf, tol=None):
    """Remove pole/zero pairs from a transfer function
    when the two match within the tolerance tol.

    tf may be a transfer function (lti object) or a list of transfer functions
    tol is optional and defaults to the system epsilon if unset.

    return a list of tfs or a tf, depending on the input type.
    """
    # initially based on python-control
    # which is in turn based on octave minreal
    # then modified considerably

    # recursively handle multiple tfs
    if not _is_zpk(tf) and not _is_num_den(tf) and not _is_A_B_C_D(tf) \
            and (isinstance(tf, list) or isinstance(tf, tuple)):
        ret = []
        for tfi in tf:
            ret += [minreal(tfi, tol)]
        return ret

    # default accuracy
    sqrt_eps = np.sqrt(eps)

    if (hasattr(tf, 'inputs') and not tf.inputs == 1) or \
       (hasattr(tf, 'outputs') and not tf.outputs == 1):
        raise TypeError("Only SISO transfer functions can be evaluated.")
    if hasattr(tf, 'zeros') and hasattr(tf, 'poles'):
        # LTI objects have poles and zeros,
        zeros = tf.zeros
        poles = tf.poles
        if hasattr(tf, 'k'):
            k = tf.k
        elif hasattr(tf, 'gain'):
            k = tf.gain
    else:
        # k = num[0] / den[0]
        zeros, poles, k = _get_zpk(tf)

    zeros = carray(zeros)
    poles = carray(poles)
    zeros.sort()
    poles.sort()
    reducedzeros = []

    for z in zeros:
        t = tol or 1000 * max(eps, abs(z) * sqrt_eps)
        idx = np.where(abs(poles - z) < t)[0]
        if len(idx):
            # cancel this zero against one of the poles
            # remove the pole and do not add the zero to the new
            poles = np.delete(poles, idx[0])
        else:
            # no matching pole
            reducedzeros.append(z)
    if len(reducedzeros):
        newzeros = carray(reducedzeros)
        num = k * np.real(np.poly(newzeros))
    else:
        num = np.array([k])
    den = np.real(np.poly(poles))

    return lti(num, den)


def diagonal_indices(a, offset=0):
    """Return the indices to the main diagonal of a 2D array a (if offset = 0),
    or to a secondary diagonal, having the offset from the main one as specified.

    The array a does not need to be square.

    Note: the sup-diagonal is at offset +1, the sub-diagonal at offset -1.

    Returns: (xs, ys)
    """
    di, dj = np.diag_indices_from(a[:min(a.shape), :min(a.shape)])
    if offset > 0:
        di, dj = zip(*[(i, j)
                     for i, j in zip(di, dj + offset) if 0 <= j < a.shape[1]])
    elif offset == 0:
        pass
    else:
        di, dj = zip(*[(i, j)
                     for i, j in zip(di - offset, dj) if 0 <= i < a.shape[0]])
    return di, dj


def circshift(a, shift):
    """Shift an array circularly.

    The ``circshift(a, shift)`` function circularly shifts the values in the
    array ``a`` by ``shift`` elements.

    **Parameters:**

    a : ndarray
        the array to be shifted. Notice that a should have a greater or equal
        number of dimensions than ``shift`` (``shift`` being a scalar is equal
        to ``shift`` being a one-dimension array.)

    shift : int or ndarray-like of int.
        the N-th element specifies the shift amount for the N-th dimension
        of the input array ``a``.

    If an element in ``shift`` is positive, the values of ``a`` are
    shifted to higher-index rows (ie down) or to higher-index columns
    (ie to the right).

    If the element is negative, the values of ``a`` are shifted in the opposite
    directions, towards lower-index rows (ie up) or to lower-index columns
    (ie right).

    If ``shift`` is an integer, the shift happens along axis 0.

    All dimensions that do not have a corresponding shift value in ``shift``
    are left untouched (ie ``shift=(1, 0, 0)`` is equal to ``shift=(1,)``,
    with the exception that the former will trigger an ``IndexError``
    if ``a.ndim < 3``).

    **Returns:**

    The shifted array.
    """
    if np.isscalar(shift):
        shift = [shift]
    for axis, ashift in enumerate(shift):
        a = np.roll(a, ashift, axis=axis)
    return a


def save_input_form(a):
    """Save the form of `a` so that it can be restored later on

    Returns: an object representing the form of `a`, to be passed to
                 restore_input_form(a, form)
    """
    if np.isscalar(a):
        ret = 'scalar'
    elif hasattr(a, 'shape'):
        ret = a.shape
    elif type(a) == tuple:
        ret = 'tuple'
    elif type(a) == list:
        ret = 'list'
    else:
        raise TypeError("Unsupported input %s" % repr(a))
    return ret


def restore_input_form(a, form):
    """Restore the form of `a` according to `form`.

    Returns: the object `a`, in the correct `form`.

        Note: use `save_input_form(a)` to get the object `form`
    """
    if form == 'scalar':
        a = np.real_if_close(a)
        if not np.isscalar(a):
            a = a.reshape((1, ))[0]
    elif form == 'tuple':
        if not type(a) == tuple:
            a = [np.real_if_close(i).reshape((1,))[0] for i in a]
            a = tuple(a)
    elif form == 'list':
        if not type(a) == list:
            a = [np.real_if_close(i).reshape((1,))[0] for i in a]
            a = list(a)
    else:
        a = a.reshape(form)
    return a


def pretty_lti(arg):
    """Given the lti object ``arg`` return a *pretty* representation."""
    z, p, k = _get_zpk(arg)
    z = np.atleast_1d(z)
    p = np.atleast_1d(p)
    z = np.round(np.real_if_close(z), 4)
    p = np.round(np.real_if_close(p), 4)
    signs = {1:'+', -1:'-'}
    if not len(z) and not len(p):
        return "%g" % k
    ppstr = ["", "", ""]
    if k != 1:
        ppstr[1] = "%g " % k
    for i, s in zip((0, 2), (z, p)):
        rz = None
        m = 1
        sorted_singularities = cplxpair(s)
        for zindex, zi in enumerate(sorted_singularities):
            zi = np.round(np.real_if_close(zi), 4)
            if np.isreal(zi):
                if len(sorted_singularities) > zindex + 1 and \
                    sorted_singularities[zindex + 1] == zi:
                    m += 1
                    continue
                if zi == 0.:
                    ppstr[i] += "z"
                else:
                    ppstr[i] += "(z %s %g)" % (signs[np.sign(-zi)], np.abs(zi))
                if m == 1:
                    ppstr[i] += " "
                else:
                    ppstr[i] += "^%d " % m
                m = 1
            else:
                if len(sorted_singularities) > zindex + 2 and \
                    sorted_singularities[zindex + 2] == zi:
                    m += .5
                    continue
                if rz is None:
                    rz = zi
                    continue
                ppstr[i] += "(z^2 %s %gz %s %g)" % \
                            (signs[np.sign(np.real_if_close(np.round(-rz - zi, 3)))],
                             np.abs(np.real_if_close(np.round(-rz - zi, 3))),
                             signs[np.sign(np.real_if_close(np.round(rz * zi, 4)))],
                             np.abs(np.real_if_close(np.round(rz * zi, 4))))
                if m == 1:
                    ppstr[i] += " "
                else:
                    ppstr[i] += "^%d " % m
                rz = None
                m = 1
        ppstr[i] = ppstr[i][:-1] if len(ppstr[i]) else "1"
    if ppstr[2] == '1':
        return ppstr[1] + ppstr[0]
    else:
        if ppstr[0] == '1' and len(ppstr[1]) and float(ppstr[1]) != 1.:
            ppstr[0] = ppstr[1][:-1]
            ppstr[1] = ""
        space_pad_ln = len(ppstr[1])
        fraction_line = "-" * (max(len(ppstr[0]), len(ppstr[2])) + 2)
        ppstr[1] += fraction_line
        ppstr[0] = " "*space_pad_ln + ppstr[0].center(len(fraction_line))
        ppstr[2] = " "*space_pad_ln + ppstr[2].center(len(fraction_line))
    return "\n".join(ppstr)


def _get_zpk(arg, input=0):
    """Utility method to convert the input arg to a z, p, k representation.

    **Parameters:**

    arg, which may be:

    * ZPK tuple,
    * num, den tuple,
    * A, B, C, D tuple,
    * a scipy LTI object,
    * a sequence of the tuples of any of the above types.

    input : scalar
        In case the system has multiple inputs, which input is to be
        considered. Input `0` means first input, and so on.

    **Returns:**

    The sequence of ndarrays z, p and a scalar k

    **Raises:**

    TypeError, ValueError

    .. warn: support for MISO transfer functions is experimental.
    """
    z, p, k = None, None, None
    if isinstance(arg, np.ndarray):
        # ABCD matrix
        A, B, C, D = partitionABCD(arg)
        z, p, k = ss2zpk(A, B, C, D, input=input)
    elif hasattr(arg, '__class__') and arg.__class__.__name__ == 'lti':
        z, p, k = arg.zeros, arg.poles, arg.gain
    elif _is_zpk(arg):
        z, p, k = np.atleast_1d(arg[0]), np.atleast_1d(arg[1]), arg[2]
    elif _is_num_den(arg):
        sys = lti(*arg)
        z, p, k = sys.zeros, sys.poles, sys.gain
    elif _is_A_B_C_D(arg):
        z, p, k = ss2zpk(*arg, input=input)
    elif isinstance(arg, collections.Iterable):
        ri = 0
        for i in arg:
            # Note we do not check if the user has assembled a list with
            # mismatched lti representations.
            if hasattr(i, 'B'):
                iis = i.B.shape[1]
                if input < ri + iis:
                    z, p, k = ss2zpk(i.A, i.B, i.C, i.D, input=input - ri)
                    break
                else:
                    ri += iis
                    continue
            elif _is_A_B_C_D(arg):
                iis = arg[1].shape[1]
                if input < ri + iis:
                    z, p, k = ss2zpk(*arg, input=input - ri)
                    break
                else:
                    ri += iis
                    continue
            else:
                if ri == input:
                    sys = lti(*i)
                    z, p, k = sys.zeros, sys.poles, sys.gain
                    break
                else:
                    ri += 1
                    continue
                ri += 1
        if (z, p, k) == (None, None, None):
            raise ValueError("The LTI representation does not have enough" +
                             "inputs: max %d, looking for input %d" %
                             (ri - 1, input))
    else:
        raise TypeError("Unknown LTI representation: %s" % arg)
    return z, p, k


def _get_num_den(arg, input=0):
    """Utility method to convert the input arg to a (num, den) representation.

    **Parameters:**

    arg, which may be:

    * ZPK tuple,
    * num, den tuple,
    * A, B, C, D tuple,
    * a scipy LTI object,
    * a sequence of the tuples of any of the above types.

    input : scalar
        In case the system has multiple inputs, which input is to be
        considered. Input `0` means first input, and so on.

    **Returns:**

    The sequence of ndarrays num, den

    **Raises:**

    TypeError, ValueError

    .. warn: support for MISO transfer functions is experimental.
    """
    num, den = None, None
    if isinstance(arg, np.ndarray):
        # ABCD matrix
        A, B, C, D = partitionABCD(arg)
        num, den = ss2tf(A, B, C, D, input=input)
    elif hasattr(arg, '__class__') and arg.__class__.__name__ == 'lti':
        num, den = arg.num, arg.den
    elif _is_num_den(arg):
        num, den = carray(arg[0]).squeeze(), carray(arg[1]).squeeze()
    elif _is_zpk(arg):
        num, den = zpk2tf(*arg)
    elif _is_A_B_C_D(arg):
        num, den = ss2tf(*arg, input=input)
    elif isinstance(arg, collections.Iterable):
        ri = 0
        for i in arg:
            # Note we do not check if the user has assembled a list with
            # mismatched representations.
            if hasattr(i, 'B'):  # lti
                iis = i.B.shape[1]
                if input < ri + iis:
                    num, den = ss2tf(i.A, i.B, i.C, i.D, input=input - ri)
                    break
                else:
                    ri += iis
            else:
                sys = lti(*i)
                iis = sys.B.shape[1]
                if input < ri + iis:
                    num, den = ss2tf(
                        sys.A, sys.B, sys.C, sys.D, input=input - ri)
                    break
                else:
                    ri += iis

        if (num, den) == (None, None):
            raise ValueError("The LTI representation does not have enough" +
                             "inputs: max %d, looking for input %d" %
                             (ri - 1, input))
    else:
        raise TypeError("Unknown LTI representation: %s" % arg)

    if len(num.shape) > 1:
        num = num.squeeze()
    if len(den.shape) > 1:
        num = den.squeeze()

    # default accuracy: sqrt_ps
    sqrt_eps = np.sqrt(eps)
    while True:
        if abs(num[0]) < sqrt_eps:
            num = num[1:]
        else:
            break
    while True:
        if abs(den[0]) < sqrt_eps:
            den = den[1:]
        else:
            break

    return num, den


def _getABCD(arg):
    """Utility method to convert the input arg to an A, B, C, D representation.

    **Parameters:**

    arg, which may be:

    * ZPK tuple,
    * num, den tuple,
    * A, B, C, D tuple,
    * a scipy LTI object,
    * a sequence of the tuples of any of the above types.

    **Returns:**

    The sequence of ndarrays A, B, C, D

    **Raises:**

    TypeError, ValueError
    """
    if isinstance(arg, np.ndarray):
        # ABCD matrix
        A, B, C, D = partitionABCD(arg)
    elif hasattr(arg, '__class__') and arg.__class__.__name__ == 'lti':
        A, B, C, D = arg.A, arg.B, arg.C, np.atleast_2d(arg.D)
    elif _is_zpk(arg) or _is_num_den(arg) or _is_A_B_C_D(arg):
        sys = lti(*arg)
        A, B, C, D = sys.A, sys.B, sys.C, sys.D
    elif isinstance(arg, collections.Iterable):
        A, B, C, D = None, None, None, None
        for i in arg:
            # Note we do not check if the user has assembled a list with
            # mismatched lti representations.
            sys = lti(*i) if not hasattr(i, 'A') else i
            if A is None:
                A = sys.A
            elif not np.allclose(sys.A, A, atol=1e-8, rtol=1e-5):
                raise ValueError("Mismatched lti list, A matrix disagreement.")
            else:
                pass
            if B is None:
                B = sys.B
            else:
                B = np.hstack((B, sys.B))
            if C is None:
                C = sys.C
            elif not np.allclose(sys.C, C, atol=1e-8, rtol=1e-5):
                raise ValueError("Mismatched lti list, C matrix disagreement.")
            if D is None:
                D = sys.D
            else:
                D = np.hstack((D, sys.D))
    else:
        raise TypeError("Unknown LTI representation: %s" % arg)
    return A, B, C, D


def _is_zpk(arg):
    """Can the argument be safely assumed to be a zpk tuple?"""
    return isinstance(arg, collections.Iterable) and len(arg) == 3 and \
        isinstance(arg[0], collections.Iterable) and \
        isinstance(arg[1], collections.Iterable) and np.isscalar(arg[2])


def _is_num_den(arg):
    """Can the argument be safely assumed to be a num, den tuple?"""
    return isinstance(arg, collections.Iterable) and len(arg) == 2 and \
        isinstance(arg[0], collections.Iterable) and \
        isinstance(arg[1], collections.Iterable)


def _is_A_B_C_D(arg):
    """Can the argument be safely assumed to be an (A, B, C, D) tuple?"""
    return isinstance(arg, collections.Iterable) and len(arg) == 4 and \
        (isinstance(arg[0], collections.Iterable) or np.is_scalar(arg[0])) and \
        (isinstance(arg[1], collections.Iterable) or np.is_scalar(arg[1])) and \
        (isinstance(arg[2], collections.Iterable) or np.is_scalar(arg[2])) and \
        (isinstance(arg[3], collections.Iterable) or np.is_scalar(arg[3]))


def _cell_like_list(shape, init=None):
    """Returns a list of lists (possibly of lists... etc...),
    with all values initialized to `init`.

    **Parameters:**

    shape: tuple of ints,
           the dimensions of the desidered object

    init: object, optional,
          the initialization value for every element in the object

    **Returns:**

    cell, a list of lists (...)

    .. deprecated:: 0.1-3

    """
    a = []
    for _ in range(shape[0]):
        if len(shape) == 1:
            a.append(init)
        else:
            a.append(_cell_like_list(shape[1:], init))
    return a


def mround(x):
    """Round ``x`` to the nearest integers.

    This is a MATLAB-compatible round function, numpy's ``round()``
    behaves differently.

    Behaviour with a fractional part of 0.5:

    * Positive elements with a fractional part of 0.5 round up to the nearest positive integer.

    * Negative elements with a fractional part of -0.5 round down to the nearest negative integer.

    If the elements of ``x`` are complex, real and imaginary parts are
    rounded separately.

    **Example:**

       >>> mround([-1.9, -0.5, -0.2, 3.4, 4.5, 5.6, 7.0, 2.4+3.6j])
       [-2.0, -1.0, 0.0, 3.0, 5.0, 6.0, 7.0, 2.0+4.0j]

    """
    iform = save_input_form(x)
    x = carray(x)

    def _mround(z):
        """Base function to generate the ufunc round"""
        z = np.real_if_close(z)
        if np.iscomplex(z):
            return _mround(np.real(z)) + 1j*_mround(np.imag(z))
        s = 1 if z >= 0 else -1
        res = z - s*(abs(z) % 1) if abs(z) % 1 < .5 \
              else z + s*(1 - (abs(z)%1))
        return res
    _internal = np.frompyfunc(_mround, 1, 1)
    xf = np.array(_internal(x), dtype=x.dtype)
    return restore_input_form(xf, iform)


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
    from ._synthesizeNTF import synthesizeNTF
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
