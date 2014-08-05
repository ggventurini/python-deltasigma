# -*- coding: utf-8 -*-
# _partitionABCD.py
# Module providing the partitionABCD function
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

"""Module providing the partitionABCD() function
"""

def partitionABCD(ABCD, m=None):
    """Partition ABCD into A, B, C, D for an m-input state-space system.

    The :math:`ABCD` matrix is defined as:

    .. math::

        ABCD =
            \\left[
            \\begin{array}{c|c}
              A & B \\\\ \\hline
              C & D
            \\end{array}
            \\right].

    The matrices :math:`A`, :math:`B`, :math:`C` and :math:`D` define the
    evolution of a generic linear discrete time invariant system under study
    through the state-space representation:

    .. math::

        x(k + 1) = Ax(k) + Bu(k)

    .. math::

        y(k) = Cx(k) + Du(k)

    The matrices are:

    * :math:`A` is the state matrix, dimensions :math:`(n, n)`, :math:`n` being the number of states in the system,

    * :math:`B` is the input matrix, dimensions :math:`(n, m)`, :math:`m` being the number of inputs in the system,

    * :math:`C` is the output matrix, dimensions :math:`(r, n)`, :math:`r` being the number of outputs in the system,

    * :math:`D` is the feedthrough matrix, dimensions :math:`(r, m)`.

    The vectors are:

    * :math:`x(k)` the state sequence, of length :math:`n`,

    * :math:`u(k)`, the input sequence, of length :math:`m`, and

    * :math:`y(k)` is the output sequence, of length :math:`r`.

    **Parameters:**

    ABCD: ndarray
          The ABCD matrix to be partitioned

    m: int, optional
       The number of inputs in the system. It will be calculated from the ABCD matrix if not provided.

    **Returns:**

    A, B, C, D: tuple of ndarrays
                The partitioned matrices.

    .. seealso::

        :ref:`Modulator model: loop filter <loop-filter-label>` for a discussion of the ABCD matrix in the particular case of a delta-sigma modulator loop filter.

    """
    # remember the ABCD matrix is assembled like this:
    # [[A, B], 
    #  [C, D]]
    if m is None:
        n = min(ABCD.shape) - 1
        m = ABCD.shape[1] - n
    else:
        n = ABCD.shape[1] - m

    r = ABCD.shape[0] - n

    A = ABCD[:n, :n]
    B = ABCD[:n, n:n+m]
    C = ABCD[n:n+r, :n]
    D = ABCD[n:n+r, n:n+m]
    return A, B, C, D

