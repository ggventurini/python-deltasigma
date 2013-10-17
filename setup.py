#!/usr/bin/env python
from distutils.core import setup
from pydelsigma import __version__

setup(name='pydelsigma',
      version=__version__,
      description='A Python translation of Richard Schreier\'s delta sigma toolbox',
      author='Giuseppe Venturini and others',
      author_email='ggventurini+GITHUB@gmail.com',
      url='http://github.com/ggventurini/python-deltasigma/',
      packages=['pydelsigma'],
      package_data={
        'pydelsigma': ['test_data/*.mat', 'test_data/*.txt']
      }
     )

print """
+---------------------------------------------------+
The following dependencies are needed for pydelsigma
to work on your system:
- numpy: http://numpy.scipy.org/
- scipy: http://www.scipy.org/
- slycot: https://github.com/repagh/Slycot
- python-control: 
  http://sourceforge.net/projects/python-control/
- matplotlib: http://matplotlib.sourceforge.net/

Recommended:
- setuptools: https://pypi.python.org/pypi/setuptools
- nose: https://pypi.python.org/pypi/nose/
- shinx: http://sphinx-doc.org/
- patience and stubbornness.

Most of which are available through PyPi.
See Install.md for more.
+---------------------------------------------------+
"""
