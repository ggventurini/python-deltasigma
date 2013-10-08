import numpy as np
from pydelsigma import evalTF

def nabsH(w, H):
	""" nabsH(w, H) computes the negative of the absolute value of H 
	at the specified frequency on the unit circle. 

	This function is used by infnorm.py.
	"""
	z = np.exp(1j*w)
	return -np.abs(evalTF(H, z))
	
def test_nabsH():
	from control.matlab import tf
	H = tf([1, 2], [2, 0, .25], 1)
	N = 129
	w = np.linspace(0, 2*np.pi, num=N, endpoint=True)
	z = np.exp(1j*w)
	r1 = -np.abs(evalTF(H, z))
	r2 = nabsH(w, H)
	assert np.allclose(r1, r2, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_nabsH()
