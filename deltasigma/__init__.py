# -*- coding: utf-8 -*-
# __init__.py
# python-deltasigma module init file
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

"""

:Author: Giuseppe Venturini
:Release: |release|
:Date: |today|

**Homepage:** http://www.python-deltasigma.io

**Documentation:** http://docs.python-deltasigma.io

**Repository:** https://github.com/ggventurini/python-deltasigma

**Bug tracker:** https://github.com/ggventurini/python-deltasigma/issues

Introduction
------------

A port of the **MATLAB Delta Sigma Toolbox** based on free software and
very little sleep.

**Python-deltasigma** is a Python package to *synthesize, simulate,
scale and map to implementable structures* **delta sigma modulators**.

It aims to provide **a 1:1 Python port** of Richard Schreier's
***excellent*** `MATLAB Delta Sigma
Toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__,
the *de facto* standard tool for high-level delta sigma simulation, upon
which it is very heavily based.\ |githalytics.com alpha|

Status
------

|Build Status| |Coverage Status| |PyPi version| |PyPi downloads| |BSD
2 clause license| |DOI BADGE|

This project is a *work in progress*, not all functionality has been
ported, yet. The next figure shows the relationship between the main functions
and the avaliable functionality at a glance.

.. image:: ../doc/_static/functionality.png

All the basic features are available since v. 0.1. A detailed changelog may be
found in `CHANGES.rst <https://github.com/ggventurini/python-deltasigma/blob/master/CHANGES.rst>`__.

Detailed information split by file and function status may be found in
`files.csv <https://github.com/ggventurini/python-deltasigma/blob/master/files.csv>`__.

The further functionality is expected to be ported and available in future releases
according to `the ROADMAP <https://github.com/ggventurini/python-deltasigma/blob/master/ROADMAP.md>`__.

Examples
--------

To see the currently implemented functionality in action, take a look
at the following ipython notebooks:

-  `dsdemo1 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo1.ipynb>`__,
   notebook port of the interactive ``dsdemo1.m``.
-  `dsdemo2 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo2.ipynb>`__,
   notebook port of the interactive ``dsdemo2.m``.
-  `dsdemo3 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo3.ipynb>`__,
   notebook port of the interactive ``dsdemo3.m``.
-  `dsdemo4 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsdemo4.ipynb>`__,
   notebook port of ``dsdemo4.m``. `Audio
   file <https://raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/sax.wav.b64>`__, right click to download.
-  `dsexample1 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsexample1.ipynb>`__, python
   version of ``dsexample1.m``.
-  `dsexample2 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsexample2.ipynb>`__, python
   version of ``dsexample2.m``.
-  `dsexample3 <http://nbviewer.ipython.org/urls/raw.githubusercontent.com/ggventurini/python-deltasigma/master/examples/dsexample3.ipynb>`__, python
   version of ``dsexample3.m``.

They are also a good means for getting started quickly.

If you have some examples you would like to share, `send me a
mail <http://tinymailto.com/5310>`__, and I will add them to the above
list.

Install
-------

Supported platforms
~~~~~~~~~~~~~~~~~~~

python-deltasigma runs on every platform and arch. supported by its
dependencies:

-  *Platforms*: Linux, Mac OS X, Windows.

-  *Archs*: x86, x86\_64 and armf (arm with floating point unit).

Dependencies
~~~~~~~~~~~~

Using ``python-deltasigma`` requires `Python 2 or 3 <http://www.python.org/>`__,
at your choice, `numpy <http://www.numpy.org/>`__,
`scipy <http://www.scipy.org>`__ (>= 0.11.0) and
`matplotlib <http://www.matplotlib.org>`__.

They are packaged by virtually all the major Linux distributions.

On a Debian Linux system, you may install them issuing:

::

    aptitude install python python-numpy python-scipy python-matplotlib

Refer to your system documentation for more information.

On Windows, I hear good things about:

-  `Enthought Canopy <https://www.enthought.com/store/>`__, a Python
   distribution that carries both free and commercial versions, and

-  `Anaconda <https://store.continuum.io/cshop/anaconda/>`__, which
   offers its full version for free.

I do not run Windows, so I can't really provide more info (sorry),
except that people tell me they manage to have a working setup.

*Mac OS X* is also supported by `Enthought
Canopy <https://www.enthought.com/store/>`__ and
`Anaconda <https://store.continuum.io/cshop/anaconda/>`__, which likely
provide the easiest and fastest solution to have a scientific Python
stack up and running.

More information can be found on the `scipy install
page <http://www.scipy.org/install.html>`__ and on the `matplotlib
homepage <http://matplotlib.org/>`__.

I wrote in a different context some directions to `compile numpy and
scipy
yourself <https://github.com/ahkab/ahkab/wiki/Install:-numpy-and-scipy>`__,
which also apply here. Be warned, it can easily get complicated.

Optional dependencies
~~~~~~~~~~~~~~~~~~~~~

The required dependencies have been kept to a minimum to allow running
``python-deltasigma`` on workstations that are not managed by the user
but by a system administrator - where typically installing libraries is
not possible and software packages are disarmingly outdated.

If at all possible, installing `Cython <http://www.cython.org>`__ is
strongly recommended.

``python-deltasigma`` contains python extension to simulate delta sigma
modulators providing a near-native execution speed -- overall roughly a
70x speed-up compared to a plain Python implementation.

On Linux, installing Cython is just one: `aptitude install cython`
away.

On Mac OS X and Windows, Cython may be installed as part of one of the
frameworks above. Please notice a compiler is needed, this may require
installing XCode and its command-line utilities or gcc through homebrew,
on Mac OS X, or Mingw, on Windows.

If the BLAS headers are found on the machine, they will be used. In
case they cannot be found automatically, it is recommended to set
the environment variable ``BLAS_H`` to the BLAS headers directory.

On Mac OS X, consider linking the headers to their conventional location::

    sudo ln -s /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/cblas.h /usr/include/cblas.h

The Cython extensions were written by Sergio Callegari, please see the
``deltasigma/`` for copyright notice and more information.

Install python-deltasigma
~~~~~~~~~~~~~~~~~~~~~~~~~

Once the dependencies set up, it is possible to install the latest stable
version directly from the `Python Package Index (PYPI)
<http://pypi.python.org>`__, running::

    pip install deltasigma

The above command will also attempt to compile and install the dependencies
in case they are not found. Please notice this is not recommended and
for this to work you should already have the required C libraries in place.

If you are interested in a bleeding-edge version -- potentially less stable
-- or in contributing code (*that's awesome!*) you can head over to
`the Github repository <http://github.com/ggventurini/python-deltasigma>`__
and check out the code from there.

Then run::

    python setup.py install

The flag ``--user`` may be an interesting option to install the package for
the current user only, and it doesn't require root privileges.

Extras for developers
~~~~~~~~~~~~~~~~~~~~~

The following may be installed at a later stage and are typically only
necessary for developers.

Building the documentation requires the
`sphinx <http://sphinx-doc.org/>`__ package. It is an optional step,
as the `the latest documentation is available
online <http://python-deltasigma.readthedocs.org/en/latest/>`__, without
need for you to build it.

If you plan to modify the code, `python-deltasigma` comes with a complete
unit tests suite, which is run against every commit and that any addition
should pass both for Python 2 and 3.
To run it locally, `setuptools <https://pypi.python.org/pypi/setuptools>`__
is needed, as it is used to access the reference function outputs.

Running the test suite may be conveniently automated installing
`nose <https://pypi.python.org/pypi/nose/>`__, and then issuing::

    nosetests -v deltasigma

from the repository root.

Useful resources
----------------

The original MATLAB Toolbox provides in-depth documentation, which is
very useful to understand what the toolbox is capable of. See
`DSToolbox.pdf <https://github.com/ggventurini/python-deltasigma/blob/master/delsig/DSToolbox.pdf?raw=true>`__
and
`OnePageStory.pdf <https://github.com/ggventurini/python-deltasigma/blob/master/delsig/OnePageStory.pdf?raw=true>`__
(*PDF warning*).

The book:

    Richard Schreier, Gabor C. Temes, *Understanding Delta-Sigma Data
    Converters*, ISBN: 978-0-471-46585-0, November 2004, Wiley-IEEE Press

is probably *the most authoritative resource on the topic*. Chapter 8-9
show how to use the MATLAB toolkit and the observations apply also to
this Python port. Links on
`amazon <http://www.amazon.com/Understanding-Delta-Sigma-Converters-Richard-Schreier/dp/0471465852>`__,
on `the Wiley-IEEE
press <http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471465852,miniSiteCd-IEEE2.html>`__.

*I am not affiliated with neither the sellers nor the authors.*

How to contribute
-----------------

Pull requests are welcome!
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to port some code, fix a bug, add a cool example or implement new
functionality (in this case it may be a good idea to get in touch early on),
*that's awesome!*

There are only a few *guidelines*, which can be overridden every time it
is reasonable to do so:

-  Please try to follow ``PEP8``.

-  Try to keep the functions signature identical. Parameters with
   ``NaN`` default values have their default value replaced with
   ``None``.

-  If a function has a varible number of return values, its Python port
   should implement the maximum number of return values.

No commit should ever fail the test suite.

Reporting bugs
~~~~~~~~~~~~~~

What bugs, *there are no bugs!*

Jokes aside, please report all bugs on `on the Github issue tracker <https://github.com/ggventurini/python-deltasigma/issues>`__.

Support python-deltasigma
~~~~~~~~~~~~~~~~~~~~~~~~~

*I do not want your money.* I develop this software because I enjoy it
and because I use it myself.

If you wish to support the development of ``python-deltasigma`` or you
find the package useful or you otherwise wish to contribute monetarily,
***please donate to cancer research instead:***

-  `Association for International Cancer Research
   (eng) <http://www.aicr.org.uk/donate.aspx>`__, or
-  `Fond. IRCCS Istituto Nazionale dei Tumori
   (it) <http://www.istitutotumori.mi.it/modules.php?name=Content&pa=showpage&pid=24>`__.

Consider `sending me a mail <http://tinymailto.com/5310>`__ afterwards,
***it makes for great motivation!***

License, copyright, rationale and credits
-----------------------------------------

Why this project was born
~~~~~~~~~~~~~~~~~~~~~~~~~

I like challenges, delta-sigma modulation and I don't have the money for
my own MATLAB license. After all, *which grad student or young researcher
has it?*

With this Python package you can simulate delta-sigma modulators for
free, on any PC.

I hope you find it useful.

Licensing and copyright notice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All original MATLAB code is Copyright (c) 2009, Richard Schreier. See
the LICENSE file for the licensing terms.

The Python code here provided is a derivative work from the above
toolkit and subject to the same license terms.

Credit goes to Richard Schreier for the original ideas, their MATLAB
implementation and the all the  diagrams found in this documentation.
Little-to-no conceptual improvements are introduced here, just code
adaptation, refactoring, rewrites and fixing of a few minor issues.

This package contains some source code from ``pydsm``, also based on the
same MATLAB toolbox. The ``pydsm`` package is copyright (c) 2012, Sergio
Callegari.

When not otherwise specified, the Python code is Copyright 2013,
Giuseppe Venturini and the python-deltasigma contributors.

MATLAB is a registered trademark of The MathWorks, Inc.

.. |githalytics.com alpha| image:: https://cruel-carlota.pagodabox.com/36f25accf60f391456efe66910bf84f8
   :target: http://githalytics.com/ggventurini/python-deltasigma
   :width: 1
.. |Build Status| image:: https://travis-ci.org/ggventurini/python-deltasigma.png?branch=master
   :target: https://travis-ci.org/ggventurini/python-deltasigma
.. |Coverage Status| image:: https://coveralls.io/repos/ggventurini/python-deltasigma/badge.png?branch=master
   :target: https://coveralls.io/r/ggventurini/python-deltasigma?branch=master
.. |PyPi version| image:: http://img.shields.io/badge/version-0.1-brightgreen.png
   :target: https://pypi.python.org/pypi/deltasigma/
.. |PyPi downloads| image::  https://pypip.in/download/deltasigma/badge.png
   :target: https://pypi.python.org/pypi/deltasigma/
.. |BSD 2 clause license| image:: http://img.shields.io/badge/license-BSD-brightgreen.png
   :target: https://raw.githubusercontent.com/ggventurini/python-deltasigma/master/LICENSE
.. |DOI BADGE| image:: https://zenodo.org/badge/doi/10.5281/zenodo.11167.png
   :target: http://dx.doi.org/10.5281/zenodo.11167


Credits
~~~~~~~

The ``python-deltasigma`` package was written by
`Giuseppe Venturini <https://github.com/ggventurini>`__, as a derivative work
of Richard Schreier's MATLAB Delta Sigma toolbox. It contains code from
``pydsm``, also based on the same MATLAB toolbox and written by Sergio
Callegari.

Contributors: Shayne Hodge

Implementation model
--------------------

The internal implementation of delta sigma modulators follows closely the one in
Richard Schreier's `MATLAB Delta Sigma
Toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__,
upon which the following documentation is very heavily based.


Modulator model
~~~~~~~~~~~~~~~

A delta-sigma modulator with a single quantizer is assumed to consist of
quantizer connected to a loop filter as shown in the diagram below.

.. image:: ../doc/_static/modulator_model.png
    :align: center
    :alt: Modulator model

.. _loop-filter-label:

The loop filter
:::::::::::::::

The loop filter is described by an :math:`ABCD` matrix. For single-quantizer
systems, the loop filter is a two-input, one-output linear system and
:math:`ABCD` is an :math:`(n+1, n+2)` matrix, partitioned into
:math:`A` :math:`(n, n)`, :math:`B` :math:`(n, 2)`, :math:`C` :math:`(1, n)`
and :math:`D` :math:`(1, 2)` sub-matrices as shown below:

.. math::

    ABCD =
        \\left[
        \\begin{array}{c|c}
          A & B \\\\ \\hline
          C & D
        \\end{array}
        \\right].

The equations for updating the state and computing the output of the loop filter are:

.. math::

    x(n + 1) = Ax(n) + B
        \\left[
        \\begin{array}{c}
          u(n) \\\\
          v(n)
        \\end{array}
        \\right]

.. math::

    y(n) = Cx(n) + D
        \\left[
        \\begin{array}{c}
          u(n) \\\\
          v(n)
        \\end{array}
        \\right].

Where :math:`u(n)` is the input sequence and :math:`v(n)` is the modulator output sequence.

This formulation is sufficiently general to encompass all single-quantizer modulators which
employ linear loop filters. The toolbox currently supports translation to/from an ABCD descrip-
tion and coefficients for the following topologies:

 * CIFB : Cascade-of-integrators, feedback form.
 * CIFF : Cascade-of-integrators, feedforward form.
 * CRFB : Cascade-of-resonators, feedback form.
 * CRFF : Cascade-of-resonators, feedforward form.
 * CRFBD : Cascade-of-resonators, feedback form, delaying quantizer.
 * CRFFD : Cascade-of-resonators, feedforward form, delaying quantizer
 * PFF : Parallel feed-forward
 * Stratos : A CIFF-like structure with non-delaying resonator feedbacks [*]_

.. [*] Contributed to the MATLAB delta sigma toolbox in 2007 by Jeff Gealow.

See :ref:`topologies-diagrams` for a block-level view of the different modulator structures.

Multi-input and multi-quantizer systems can also be described with an
ABCD matrix and the previous equation
will still apply. For an :math:`n_i`-input, :math:`n_o`-output modulator,
the dimensions of the sub-matrices are
:math:`A`: :math:`(n, n)`, :math:`B`: :math:`(n, n_i + n_o)`,
:math:`C`: :math:`(n_o, n)` and :math:`D`: :math:`(n_o, n_i+n_o)`.

Quantizer model
:::::::::::::::

The quantizer is ideal, producing integer outputs centered about zero. Quantizers with an even
number of levels are of the mid-rise type and produce outputs which are odd integers. Quantizers
with an odd number of levels are of the mid-tread type and produce outputs which are even inte-
gers.

.. image:: ../doc/_static/quantizer_model.png
    :align: center
    :alt: Quantizer model

.. seealso:: :func:`bquantize`, :func:`bunquantize`

.. _topologies-diagrams:

Topologies diagrams
~~~~~~~~~~~~~~~~~~~

All the following topology diagrams are reproduced from
`DSToolbox.pdf <https://github.com/ggventurini/python-deltasigma/blob/master/delsig/DSToolbox.pdf?raw=true>`__
in the `MATLAB Delta Sigma Toolbox <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__,
written by Richard Schreier. All credits belong to the original author.

CIFB
::::

.. image:: ../doc/_static/CIFB.png
    :align: center
    :alt: CIFB topology

CIFF
::::

.. image:: ../doc/_static/CIFF.png
    :align: center
    :alt: CIFF topology

CRFB
::::

.. image:: ../doc/_static/CRFB.png
    :align: center
    :alt: CRFB topology

CRFF
::::

.. image:: ../doc/_static/CRFF.png
    :align: center
    :alt: CRFF topology

CRFBD
:::::

.. image:: ../doc/_static/CRFBD.png
    :align: center
    :alt: CRFBD topology

CRFFD
:::::

.. image:: ../doc/_static/CRFFD.png
    :align: center
    :alt: CRFFD topology

.. _discrete-time-to-continuous-time-mapping:

Discrete time to continuous time mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The approach presented here to design a CT delta-sigma modulator starts with the
synthesis of the noise transfer function.

Once a suitable NTF has been identified, we need to realize the transfer function with a
continuous time loop filter.

First, a loop filter topology is selected among the Feed-Forward (FF) and Feedback (FB)
structures. The feed-forward or feedback paths -- depending on the topology -- will be
characterized by an unknown proportionality factor :math:`k_i`, for each of the
:math:`i \in \\{0 \dots order\\}` branches.

.. image:: ../doc/_static/DS_equivalence_DT_CT.png
    :align: center
    :alt: DT-CT mapping model

The objective is to determine the gain factors :math:`{k_i}` to construct a CT loop
filter such that its sampled pulse response is equal to the impulse response of a DT
prototype's loop filter.

We consider here three approaches:

* equating the loop filter transfer functions,
* matching the loop filter impulse responses, implemented in :func:`realizeNTF_ct`
  as method ``'LOOP'``,
* matching the DT NTF and CT filter pulse responses, implemented in
  :func:`realizeNTF_ct` as method ``'NTF'``.

1: equate the loop filter transfer functions
::::::::::::::::::::::::::::::::::::::::::::

The DT loop filter transfer function can be found from:

.. math::

    L_{1,DT}(z) = 1 - \\frac{1}{NTF_{DT}(z)}

The CT loop filter is readily known, since it has been selected by the user,
but it needs to be converted to an equivalent DT transformation. This operation
is performed through the impulse invariance transformation, here denoted as
:math:`\\mathscr{I\\!I\\!T}` [R2]_.

Solving the equation:

.. math::

    L_{1, DT}(z) = \\mathscr{I\\!I\\!T}\\!\\!\\left(L_{1, CT}\\right)(z)

will allow determining the exact values of the coefficients :math:`k_i` [R1]_.

This approach has limited use in the real case:

* Developing an analytical CT model and applying the impulse invariance
  transformation may entail significant difficulties in presence of
  non-idealities.
* The equation above will not have a solution when the integrators are non-ideal
  since the poles of the transfer functions on RHS and LHS are different.

2: match the loop filter impulse responses
::::::::::::::::::::::::::::::::::::::::::

Although it may be non-trivial to reach an analytical expression for
:math:`\\mathscr{I\\!I\\!T}\\!\\!\\left(L_{1, CT}\\right)(z)`, it is still
possible to evaluate through numerical simulations (or in a some case
analytically), N samples of the pulse response of each path
:math:`\\{l_i[n]\\}`, obtained setting all the gain coefficients to one
or :math:`k_i = 1, \\; \\forall i \\in \\{0 \\dots order\\}`.

Evaluating the impulse response of the prototype DT loop filter, denoted in the
following as :math:`l[n]`, it is possible to write the equation:

.. math::

    [\\ l_0[n]\\; l_1[n]\\; \\dots \\; l_{order}[n]\\ ]\ K = l[n]

Where we define the vector :math:`K` as:

.. math::

    K = [\\ k_0 \\; k_1 \\; \dots \\; k_{order}\\ ]^T

In the ideal case, provided that the impulse responses have been evaluated for
a sufficiently high number of points :math:`N` (:math:`N > order`), the
equation has a exact solution, indepently of :math:`N` [R2]_.

In presence of non-idealities, it is possible to use Least Squares fitting to
find the optimum :math:`\{k_i\}`.

As discussed in [R3]_, this method is particularly sensitive to the value of
:math:`N` in the non-ideal cases.

3: match the DT NTF and CT filter pulse responses
:::::::::::::::::::::::::::::::::::::::::::::::::

A more robust approach is the following:

Remember the goal of DT to CT mapping is to have:

.. math::

    NTF_{DT} = NTF_{eq, CT}

And by definition:

.. math::

    NTF = \\frac{1}{1+L_1(z)}

Combining the two we can write the equation:

.. math::

    NTF_{DT} = \\frac{1}{1+L_{1,eq,CT}(z)}

To avoid having to apply the impulse invariance transformation, we can
rewrite the above in the time domain, getting:

.. math::

    \\sum_i k_i \\left(h[n] * l_i[n]\\right) = \\delta[n] - h[n] \\qquad \\forall i

Where :math:`h[n]` is the impulse response of the DT NTF.

The above can be rewritten as:

.. math::

    [\\ h[n]*l_0[n] \\; \\dots \\; h[n]*l_{order}[n]\\ ] K = \\delta[n] - h[n]

And solved exactly, in the ideal case, or in the least squares sense,
in presence of non-idealities [R3]_.

.. [R1] P. Benabes, M. Keramat, and R. Kielbasa, "A methodology for designing
        continuous-time sigma-delta modulators," in Proc. Eur. Conf. Design
        Test, Washington, DC, 1997, pp. 46-50

.. [R2] J. Cherry and W. Snelgrove, "Excess loop delay in continuous-time
        delta-sigma modulators," IEEE Trans. Circuits Syst. II,
        Analog Digit. Signal Process., vol. 46, no. 4, pp. 376-389, Apr. 1999

.. [R3] Shanthi Pavan, "Systematic Design Centering of Continuous Time
        Oversampling Converters",  IEEE Trans. on Circuits Syst. II Express
        Briefs, Volume.57, Issue 3, pp.158, 2010

Package contents
----------------

Key Functions
~~~~~~~~~~~~~

.. autosummary::
    :nosignatures:

    synthesizeNTF
    clans
    synthesizeChebyshevNTF
    simulateDSM
    simulateSNR
    realizeNTF
    stuffABCD
    mapABCD
    scaleABCD
    calculateTF
    realizeNTF_ct
    mapCtoD
    evalTFP

Functions for quadrature systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Notice the current version of ``python-deltasigma`` cannot synthesize
quadrature modulators. This feature is expected in v. 0.2.

Nontheless, the following functions are provided since if you already have
a synthesized modulator (or know its ABCD matrix), they allow you to simulate
the modulator with the standard tools for real ABCDs topologies.

.. autosummary::
    :nosignatures:

    mapQtoR
    mapRtoQ

Other selected functions
~~~~~~~~~~~~~~~~~~~~~~~~

The following are auxiliary functions that complement the key functions above.

.. autosummary::
    :nosignatures:

    mod1
    mod2
    calculateSNR
    predictSNR
    partitionABCD
    infnorm
    impL1
    l1norm
    pulse
    rmsGain

Utility functions for simulation of delta-sigma modulators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Functions for low-level handling of delta-sigma modulator representations,
their evaluation and filtering.

.. autosummary::
    :nosignatures:

    bquantize
    bunquantize
    cancelPZ
    circshift
    delay
    ds_f1f2
    ds_freq
    ds_hann
    ds_optzeros
    ds_quantize
    ds_synNTFobj1
    dsclansNTF
    evalMixedTF
    evalRPoly
    evalTF
    nabsH
    peakSNR
    sinc_decimate
    zinc

General utilities for data processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following are generic functions, useful for misc. tasks, like
manipulating data, conversions or padding, for example.
They provide speciality functions which are not otherwise available
in the usual scientific Python stack.

.. autosummary::
    :nosignatures:

    db
    dbm
    dbp
    dbv
    undbm
    undbp
    undbv
    rms
    padb
    padl
    padr
    padt
    cplxpair
    mfloor
    mround
    rat
    gcd
    lcm

Plotting and data display utilitites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Graphic functions:

.. autosummary::
    :nosignatures:

    plotPZ
    plotSpectrum
    figureMagic
    DocumentNTF
    PlotExampleSpectrum
    axisLabels
    bilogplot
    lollipop
    changeFig

Textual and non-graphic display-related functions

.. autosummary::
    :nosignatures:

    circ_smooth
    bplogsmooth
    logsmooth
    pretty_lti
    SIunits

All functions in alphabetical order
-----------------------------------

"""

__author__ = "Giuseppe Venturini and the python-deltasigma contributors"
__copyright__ = "Copyright 2013, Giuseppe Venturini"
__credits__ = ["Giuseppe Venturini"]
__license__ = "BSD 2-Clause License"
__version__ = '0.1-6'
__maintainer__ = "Giuseppe Venturini"
__email__ = "ggventurini+github@gmail.com"
__status__ = "Stable"

# Package testing can be done remotely, without display. This would make
# matplotlib fail (and consequently, the test itself).
# We check for $DISPLAY, but this makes us probably lose in portability,
# does Windows have this environment variable defined?
# Then again, who runs the test suit on a head-less windows machine...
# ... for the time being the following should be OK. If in the future that
# feature is needed by somebody, we can switch to
# if not os.system('python -c "import matplotlib.pyplot as plt;plt.figure()"')
import matplotlib
import os
if not ('DISPLAY' in os.environ or os.environ.get('READTHEDOCS', None)):
    matplotlib.use('Agg')

from ._DocumentNTF import DocumentNTF
from ._PlotExampleSpectrum import PlotExampleSpectrum
from ._SIunits import SIunits
from ._axisLabels import axisLabels
from ._bilogplot import bilogplot
from ._bplogsmooth import bplogsmooth
from ._bquantize import bquantize
from ._bunquantize import bunquantize
from ._calculateSNR import calculateSNR
from ._calculateTF import calculateTF
from ._cancelPZ import cancelPZ
from ._changeFig import changeFig
from ._circ_smooth import circ_smooth
from ._clans import clans
from ._constants import eps
from ._db import db
from ._dbm import dbm
from ._dbp import dbp
from ._dbv import dbv
from ._delay import delay
from ._ds_f1f2 import ds_f1f2
from ._ds_freq import ds_freq
from ._ds_hann import ds_hann
from ._ds_optzeros import ds_optzeros
from ._ds_quantize import ds_quantize
from ._ds_synNTFobj1 import ds_synNTFobj1
from ._dsclansNTF import dsclansNTF
from ._evalMixedTF import evalMixedTF
from ._evalRPoly import evalRPoly
from ._evalF0 import evalF0
from ._evalF1 import evalF1
from ._evalTF import evalTF
from ._evalTFP import evalTFP
from ._figureMagic import figureMagic
from ._frespF1 import frespF1
from ._impL1 import impL1
from ._infnorm import infnorm
from ._l1norm import l1norm
from ._logsmooth import logsmooth
from ._lollipop import lollipop
from ._mapABCD import mapABCD
from ._mapCtoD import mapCtoD
from ._mapQtoR import mapQtoR
from ._mapRtoQ import mapRtoQ
from ._mod1 import mod1
from ._mod2 import mod2
from ._nabsH import nabsH
from ._padb import padb
from ._padl import padl
from ._padr import padr
from ._padt import padt
from ._partitionABCD import partitionABCD
from ._peakSNR import peakSNR
from ._plotPZ import plotPZ
from ._plotSpectrum import plotSpectrum
from ._predictSNR import predictSNR
from ._pulse import pulse
from ._realizeNTF import realizeNTF
from ._realizeNTF_ct import realizeNTF_ct
from ._rms import rms
from ._rmsGain import rmsGain
from ._scaleABCD import scaleABCD
from ._simulateDSM import simulateDSM, simulation_backends
from ._simulateSNR import simulateSNR
from ._sinc_decimate import sinc_decimate
from ._stuffABCD import stuffABCD
from ._synthesizeChebyshevNTF import synthesizeChebyshevNTF
from ._synthesizeNTF import synthesizeNTF
from ._thermometer import thermometer
from ._undbm import undbm
from ._undbp import undbp
from ._undbv import undbv
from ._utils import circshift, cplxpair, mfloor, mround, pretty_lti, rat, gcd, lcm
from ._zinc import zinc
