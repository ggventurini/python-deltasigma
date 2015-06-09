CHANGES
~~~~~~~

Version 0.2 series
------------------

The 0.2 series features support for real and quadrature baseband and passband
modulator topologies.

**0.2.0**: Add support for quadrature modulators.
 * Add ``simulateQDSM()``,
 * Add ``synthesizeQNTF()``,
 * Add ``realizeQNTF()``,
 * Add ``simulateQSNR()``,
 * Add ``calculateQTF()``.
 * Several functions have been extened to support quadrature modulator.

Version 0.1 series
------------------

**0.1-10**: Bugfix and additional options in ``changeFig()``.
 * BUGFIX: ship Cython sources in the package.
 * BUGFIX: only modify explicitely set options in ``changeFig()``.
 * Add support for ``xfticks`` and ``yfticks`` options in ``changeFig()``.
 * Add support for BW conversion of plots in ``changeFig()``.

**0.1-9**: Add support for modulators with multiple outputs, allowing simulating MASH cascade DSMs.
 * Add support in ``partitionABCD()`` for specifying the number of outputs.
 * Add support for multiple quantizers in ``calculateTF()``.
 * BUGFIX: Fix simulation of DSMs with multiple quantizers.
 * BUGFIX: ``cancelPZ()`` was not testing the first root.
 * ``pretty_lti()`` now returns 0 if k is 0 after rounding.
 * ``plotPZ()`` doesn't list coincident real roots as complex conjugate
   roots with imag(root) = +/-0.
 * DOC: add example of MASH cascade.

**0.1-8**: Accept both tuples and lists as NTFs for simulation.
 * Previously passing a tuple for the NTF resulted in an error. Fixed.
 * Doc fixes in ``mapCtoD()``.

**0.1-7**: Quantization fix in the Cython backends. More tests.
 * A bug was found in the function responsible for quantizing the loop
   filter output in ``simulateDSM()``, only the Cython implementations are
   affected: **all users are strongly recommended to upgrade**.
 * Add more test for ``simulateDSM()``.
 * Check for the filter and data lengths in ``sinc_decimate()``.

**0.1-6**: ``sinc_decimate()`` fix, NTF matching method in ``realizeNTF_ct()``
 * An off-by-1 indexing bug was found in ``sinc_decimate()``, **all users are
   strongly recommended to update**.
 * Add NTF matching method to ``realizeNTF_ct()``.

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
