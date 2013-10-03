import numpy as np

def undbm(p, z=50):
	""" v = undbm(p, z=50) = sqrt(z*10^(p/10-3))
	RMS voltage equivalent to a power p in dBm
	z is the normalization resistance (defaults to 50ohm)
	"""
	return np.sqrt(z*10.**(p/10.-3))

def test():
	assert np.allclose([undbm(53.015)], [100.054125892], rtol=1e-05, atol=1e-08)
	assert np.allclose([undbm(3, 100)], [0.44668359215], rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()
