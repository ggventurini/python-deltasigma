# utils.py
# Miscellaneous functions and wrappers for MATLAB functions 
# that do not find a direct replacement in numpy

import numpy as np
import fractions
from fractions import Fraction as Fr

def rat(x, tol):
	"""num, den = rat(x, tol)
	where num/den == x to the specified tolerance tol
	Note: num, den are of type 'int'
	"""
	return Fr(float(x)).limit_denominator(int(1/float(tol))).numerator, \
                Fr(float(x)).limit_denominator(int(1/float(tol))).denominator

gcd = fractions.gcd

lcm = lambda a, b: int(a*b / float(gcd(a, b)))

def test_rat():
	import numpy.random as rnd
	for i in range(10):
		n, d = rnd.randint(1, 5000), rnd.randint(1, 5000)
		fr = float(n)/float(d)
		nt, dt = rat(fr, tol=200e-6)
		assert np.allclose((n/float(d) -  nt/float(dt),), (0.,), atol=200e-6, rtol=1e-12)

def test_gcd_lcm():
	a, b = 36, 31721
	tlcm, tgcd = 1141956, 1
	assert lcm(a, b) == tlcm
	assert gcd(a, b) == tgcd

if __name__ == '__main__':
	test_rat()
	test_gcd_lcm()
