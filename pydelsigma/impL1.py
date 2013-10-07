from pydelsigma import padr
import numpy as np
import scipy
from scipy.signal import convolve as conv
from control.matlab import impulse, tf

# man this is ewwww
padr = padr.padr

def impL1(arg1, n=10):
	""" y = impL1(ntf, n=10) 
	Compute the impulse response from the comparator
	output to the comparator input for the given NTF.
	n is the (optional) number of points (10).
	 
	This function is useful when verifying the realization
	of a NTF with a specified topology.
	"""
	if not hasattr(arg1, 'pole'):
		raise ValueError, 'No poles field in the NTF.'
	if not hasattr(arg1, 'zero'):
		raise ValueError, 'No zeros field in the NTF.'

	z = arg1.zero()
	p = arg1.pole()

	lf_den = padr(arg1.num[0][0], len(p)+1)
	lf_num = lf_den - arg1.den[0][0]
	ts = np.arange(n)
	all_lf = np.concatenate((lf_num, lf_den), axis=1)
	if not np.allclose(np.imag(all_lf), np.zeros(all_lf.shape), atol=1e-9):
		# Complex loop filter
		lfr_den = np.real(conv(lf_den, np.conj(lf_den)))
		lfr_num = conv(lf_num, np.conj(lf_den))
		lf_i = tf(np.real(lfr_num).tolist()[0], lfr_den.tolist()[0], 1)
		lf_q = tf(np.imag(lfr_num).tolist()[0], lfr_den.tolist()[0], 1)
		y = impulse(lf_i, T=ts) + 1j * impulse(lf_q, T=ts)
	else:
		y = impulse(tf(lf_num.tolist()[0], lf_den.tolist()[0], 1), T=ts)
	return y

def test_impL1():
	"""Test function for impL1()"""
	sys1 = tf([1], [1, 2, 1])
	r1, t1 = impL1(sys1, n=10)
	r2, t2 = np.array([[ -2.00000000e+00, -1.00000000e+00, 2.41009246e-16, 1.95887461e-17, 
		    -8.67469990e-33, -3.83718972e-34, 2.47373167e-49, 7.51657350e-51, 
		    -6.36281335e-66, -1.47240249e-67]]), np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
	assert np.allclose(r1, r2, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_impL1()
