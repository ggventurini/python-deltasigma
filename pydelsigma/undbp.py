import numpy as np

def undbp(x):
	""" y = undbp(x)	
	Convert x from dB to a power"""
	return 10.**(x/10.)
	
def test():
	assert np.allclose([undbp(53.05)], [201836.636368], rtol=1e-05, atol=1e-08)
	assert np.allclose([undbp(3)], [1.99526231497], rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()
