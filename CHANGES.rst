CHANGES
~~~~~~~

Version 0.1 series
------------------

The 0.1 series features support for (real) baseband and passband modulator
topologies.

**0.1-5**: CRFFD support, separate tests, less verbosity and DOC fixes.
 * Add CRFFD support (see ``realizeNTF``, ``mapABCD`` and ``stuffABCD``).
 * Move all tests to a dedicated location (tests/).
 * Ensure float64 is the data representation when simulating DSM.
 * Add the ``simulations_backends`` variable and its doc.
 * Cython: disable cblas extension on Win. Reduce verbosity.
 * Multiple minor fixes to ensure scalars are never returned in place of
   arrays.

Many thanks to Shayne Hodge for reporting issues with ``deltasigma`` on
Windows and several patches to the test suite.

**0.1-4**: Cython implementation of ``simulateDSM()``, PEP8 and DOC fixes.
 * ``deltasigma/_simulateDSM_cblas.pyx`` and
   ``deltasigma/_simulateDSM_scipy_blas.pyx``, Cython implementation from
   ``pydsm`` of ``simulateDSM()``, available if Cython is, providing a 70x
   speed-up of DSM simulations.
 * More documentation improvements and PEP8-related fixes.

**0.1-3**: Bugfixes, PEP8, more test coverage
 * ``deltasigma/_realizeNT_ct.py`` now supports ``FF`` topologies.
 * ``deltasigma/_pulse.py`` now supports MIMO systems.
 * ``pretty_lti()`` has been improved to provide the prettiest printing of LTIs
   to date.
 * Many documentation improvements and PEP8-related fixes.

**0.1-2**: Bugfixes, PEP8, DOC and most importantly a, g, b, c reshape.
 * The a, g, b, c coefficients are now 1-dimensional.
 * ``deltasigma/_stuffABCD.py``: scalar ``b`` bugfix.
 * ``deltasigma/_logsmooth.py``: fix bin width.
 * ``deltasigma/_utils.py``: add ``mround()``, round compatibly with MATLAB.
 * ``deltasigma/_utils.py``: add root multiplicity support in ``pretty_lti()``.
 * ``deltasigma/_utils.py``: bugfix in cplxpair for incoherent complex values.

**0.1-1**: Bugfix: most importantly fix ``realizeNTF_ct()``.
 * ``deltasigma/_realizeNTF_ct.py``: Fixes for multi-timing, add unit tests for FB.
 * ``deltasigma/_pulse.py``: Bugfix (reshape missing assignment), fix documentation formatting.
 * ``deltasigma/_bilogplot.py``: Fix plot. Add unit test.
 * ``deltasigma/_rmsGain.py``: Fix docstring.
 * ``deltasigma/_lollipop.py``: Use matplotlib's stem function. Enforce PEP8.
   Add support for color 'None'.

0.1: Bugfix: missing ``copy()`` in ``mapABCD()``.

0.1rc4 : Multiple bugfixes. Py3k fixes. Test coverage up to 85+%.

0.1rc3 : Fix file-not-found issue with ``setup.py``.

0.1rc2 : Fix travis and coveralls.io support.

0.1rc1 : Initial release
