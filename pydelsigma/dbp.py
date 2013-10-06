import numpy as np

def dbp(x):
	""" dbp(x) = 10*log10(x): the dB equivalent of the power x"""
	if not len(x):
		return 
	y = -np.inf*np.ones(x.shape)
	nonzero = (x != 0)
	y[nonzero] = 10.*np.log10(np.abs(x[nonzero]))
	return y

def test_dbp():
	tv = np.array([2])
	r = np.array([3.01029996])
	res = dbp(tv)
	assert np.allclose(r, res, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_dbp()
