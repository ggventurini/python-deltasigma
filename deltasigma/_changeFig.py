# -*- coding: utf-8 -*-
# _changeFig.py
# Modify the size of lines, markers and text labels in a plot.
# Supports subplots too.
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

"""Module providing the changeFig() function
"""

import matplotlib as mpl
import pylab as plt
import numpy as np

def changeFig(fontsize=9, linewidth=1, markersize=6, fig=None):
    """Quickly change several figure parameters.

    Take each axes in the figure, and for each line and text item 
    in the axes, set linewidth, markersize and font size.

    **Parameters:**

    fontsize : scalar, optional
        the font size, given in points

    linewidth : scalar, optional
        the line width, given in points

    markersize : scalar, optional
        the marker size, given in points

    fig : a matplotlib figure object, optional
        the figure to apply the modifications to, if not given it is assumed
        to be a the currently active figure.

    **Returns:**

    None.

    .. note:: 

        This function may be useful to enhance the readibility of
        figures to be used in presentations.

    """
    if fig is None:
        fig = plt.gcf()
    for ax in fig.get_axes():
        _setAxLinewidth(ax, linewidth, markersize)
        _setAxLabelsFontsize(ax, fontsize)
        _setTextFontsize(ax, fontsize)

def _setAxLinewidth(ax, linewidth, markersize, BW=False):
    """Take each Line2D in the axes, ax, and convert the line style
    Optionally also convert to BW, from http://tinyurl.com/qylqgoz
    """
    COLORMAP = {
        'b': {'marker': None, 'dash': (None,None)},
        'g': {'marker': None, 'dash': [5,5]},
        'r': {'marker': None, 'dash': [5,3,1,3]},
        'c': {'marker': None, 'dash': [1,3]},
        'm': {'marker': None, 'dash': [5,2,5,2,5,10]},
        'y': {'marker': None, 'dash': [5,3,1,2,1,10]},
        'k': {'marker': 'o', 'dash': (None,None)} #[1,2,1,10]}
        }

    for line in ax.get_lines():
        if BW:
            origColor = line.get_color()
            line.set_color('black')
            line.set_dashes(COLORMAP[origColor]['dash'])
            line.set_marker(COLORMAP[origColor]['marker'])
        if markersize is not None:
            line.set_markersize(markersize)
        if markersize is not None:
            line.set_linewidth(linewidth)

def _setAxLabelsFontsize(ax, fontsize):
    """Change the font size of the axis labels
    """
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)

def _setTextFontsize(ax, fontsize):
    """Change the font size of the text Artists
    """
    for artist in ax.get_children():
        if isinstance(artist, mpl.text.Text):
            artist.set_size(fontsize)

def test_changeFig():
    """Test function for changeFig()"""
    fig = plt.figure()
    xval = np.arange(0, 1, .01)

    ax = fig.add_subplot(211)
    ax.plot(xval,np.cos(2*np.pi*xval))
    ax = fig.add_subplot(212)
    ax.plot(xval,np.cos(8*np.pi*xval))
    ax.text(.01, -.9, "This is a TEST")
    changeFig(linewidth=3, fontsize=18)

    fig = plt.figure()
    xval = np.arange(0, 1, .01)

    ax = fig.add_subplot(211)
    ax.plot(xval,np.cos(2*np.pi*xval))
    ax.plot(xval,np.cos(3*np.pi*xval))
    ax.plot(xval,np.cos(4*np.pi*xval))
    ax.plot(xval,np.cos(5*np.pi*xval))
    ax.plot(xval,np.cos(6*np.pi*xval))
    ax.plot(xval,np.cos(7*np.pi*xval))
    ax.plot(xval,np.cos(8*np.pi*xval))

    ax = fig.add_subplot(212)
    ax.plot(xval,np.cos(2*np.pi*xval))
    ax.plot(xval,np.cos(3*np.pi*xval))
    ax.plot(xval,np.cos(4*np.pi*xval))
    ax.plot(xval,np.cos(5*np.pi*xval))
    ax.plot(xval,np.cos(6*np.pi*xval))
    ax.plot(xval,np.cos(7*np.pi*xval))
    ax.plot(xval,np.cos(8*np.pi*xval))
    ax.text(.01, -.9, "This is a TEST")

    #fig.savefig("origDemo.png")
    changeFig(linewidth=3, fontsize=18)
    _setAxLinewidth(ax, linewidth=3, markersize=1, BW=True)
    #fig.savefig("changeDemo.png")
    assert True #no errors
