#!/usr/bin/env python
from distutils.core import setup
from pydelsigma import __version__

setup(name='pydelsigma',
      version=__version__,
      description='A Python translation of Richard Schreier\'s delta sigma toolbox',
      author='Giuseppe Venturini and others',
      author_email='ggventurini+GITHUB@gmail.com',
      url='http://github.com/ggventurini/python-delsigma/',
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
- matplotlib: http://matplotlib.sourceforge.net/

The unit tests require scipy: http://www.scipy.org/

They are available through PyPi. (See Install.md)
+---------------------------------------------------+
"""
