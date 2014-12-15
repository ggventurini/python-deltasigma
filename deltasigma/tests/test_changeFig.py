# test_changeFig.py
# This module provides the tests for the CHANGEME function.
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

"""This module provides the test class for the changeFig() function.
"""

import unittest
import numpy as np
import pylab as plt
import deltasigma as ds

from deltasigma import changeFig

class TestchangeFig(unittest.TestCase):
    """Test class for changeFig()"""

    def setUp(self):
        pass

    def test_changeFig(self):
        """Test function for changeFig()"""
        fig = plt.figure()
        xval = np.arange(0, 1, .01)

        ax = fig.add_subplot(211)
        ax.plot(xval, np.cos(2*np.pi*xval))
        ax = fig.add_subplot(212)
        ax.plot(xval, np.cos(8*np.pi*xval))
        ax.text(.01, -.9, "This is a TEST")
        changeFig(linewidth=3, fontsize=18, xfticks='sci', yfticks='sci')

        fig = plt.figure()
        xval = np.arange(0, 1, .01)

        ax = fig.add_subplot(211)
        ax.plot(xval, np.cos(2*np.pi*xval))
        ax.plot(xval, np.cos(7*np.pi*xval))
        ax.plot(xval, np.cos(8*np.pi*xval))

        ax = fig.add_subplot(212)
        ax.plot(xval, np.cos(3*np.pi*xval))
        ax.plot(xval, np.cos(4*np.pi*xval))
        ax.plot(xval, np.cos(6*np.pi*xval))
        ax.text(.01, -.9, "This is a TEST")

        #fig.savefig("origDemo.png")
        changeFig(linewidth=3, fontsize=18, markersize=1, bw=True)
        changeFig(xfticks='plain', yfticks='plain')
        ds._changeFig._setAxLinewidth(ax, linewidth=3, markersize=1, BW=True)
        #fig.savefig("changeDemo.png")
        self.assertTrue(True) #no errors

