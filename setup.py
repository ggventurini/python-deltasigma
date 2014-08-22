#!/usr/bin/env python

import os
from setuptools import setup, find_packages
__version__ = "0.1-6"

def read(fname):
    try:
        with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
            return fp.read()
    except IOError:
        return ""

setup(
    name='deltasigma',
    version=__version__,
    packages=find_packages(),
    package_data={
      'deltasigma': ['tests/test_data/*.mat', 'tests/test_data/*.txt']
    },
    install_requires=['numpy', 'scipy>=0.11.0', 'matplotlib>=1.1.1'],
    zip_safe=False,
    include_package_data=True,
    author="Giuseppe Venturini and others",
    author_email="giuseppe.g.venturini@ieee.org",
    description="a Python package to synthesize, simulate, scale and map " + \
                "to implementable topologies delta sigma modulators.",
    long_description=''.join([read('pypi_description.rst'), '\n\n',
                              read('CHANGES.rst')]),
    license="BSD",
    keywords="delta sigma modulator simulator",
    url="http://github.com/ggventurini/python-deltasigma",
    test_suite = "deltasigma.tests",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Natural Language :: English",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3"]
)

