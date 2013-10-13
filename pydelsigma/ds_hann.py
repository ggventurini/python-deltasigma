import numpy as np

def ds_hann(n):
	""" function w = ds_hann(n)
	A Hann window of length n. Does not smear tones located exactly in a bin.
	"""
	x = np.arange(n, dtype='float_')
	return .5*(1 - np.cos(2*np.pi*x/n))

def test_ds_hann():
	res = np.array([ 0.        ,  0.02148628,  0.06768441,  0.0954915 ,  0.06533781,
	                -0.03015369, -0.1545085 , -0.24133259, -0.22851372, -0.0954915 ])
	assert np.allclose(res, np.hanning(10) - ds_hann(10), atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_ds_hann()
