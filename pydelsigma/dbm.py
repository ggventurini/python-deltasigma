from __future__ import division
import numpy as np

def dbm(v, R=50):
	""" dbm(v, R=50) = 10*log10(v^2/50*1000)  
	The equivalent in dBm of an rms voltage v
	"""
	if not len(v):
		return
	
	y = -np.Inf*np.ones(np.size(v))
	nonzero = (v != 0)
	y[nonzero] = 10.*np.log10(np.abs(v[nonzero]**2.)/R) + 30
	return y

def test_dbm():
	v = np.arange(10)*1e-3
	r = [-np.inf, -46.98970004, -40.96910013, -37.44727495, -34.94850022,
	     -33.01029996, -31.42667504, -30.08773924, -28.9279003, -27.90484985]
	assert np.allclose(dbm(v), r, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_dbm()
