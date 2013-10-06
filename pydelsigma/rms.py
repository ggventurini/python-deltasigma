import numpy as np
import numpy.linalg as la

def rms(x, no_dc=False):
	""" y = rms(x, no_dc=False)
	Calculates the Root Mean Square of the input vector.
	If no_dc is True (or non-zero), the DC value gets subtracted.
	"""
	if no_dc:
	    x = x - np.mean(x)
	return la.norm(x)/np.sqrt(len(x))
	
def test():
	tv = np.arange(100)
	res1 = np.sqrt(np.sum(tv**2.)/float(tv.shape[0]))
	res2 = np.sqrt((np.sum((tv - tv.mean())**2.))/tv.shape[0])
	#print res1, res2, rms(tv), rms(tv, True)
	assert np.allclose(rms(tv), res1, rtol=1e-05, atol=1e-08)
	assert np.allclose(rms(tv, no_dc=True), res2, rtol=1e-05, atol=1e-08)
	
if __name__ == '__main__':
	test()

