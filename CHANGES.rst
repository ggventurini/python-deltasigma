CHANGES
~~~~~~~

Version 0.1 series
------------------

The 0.1 series features support for (real) baseband and passband modulator
topologies.

**0.1-1**: Bugfix: most importantly fix ``realizeNTF_ct()``.

 * ``deltasigma/_realizeNTF_ct.py``: Fixes for multi-timing, add unit tests for FB.
 * ``deltasigma/_pulse.py``: Bugfix (reshape missing assignment), fix documentation formatting.
 * ``deltasigma/_bilogplot.py``: Fix plot. Add unit test.
 * ``deltasigma/_rmsGain.py``: Fix docstring.
 * ``deltasigma/_lollipop.py``: Use matplotlib's stem function. Enforce PEP8. Add support for color 'None'.

0.1: Bugfix: missing ``copy()`` in ``mapABCD()``.

0.1rc4 : Multiple bugfixes. Py3k fixes. Test coverage up to 85+%.

0.1rc3 : Fix file-not-found issue with ``setup.py``.

0.1rc2 : Fix travis and coveralls.io support.

0.1rc1 : Initial release
