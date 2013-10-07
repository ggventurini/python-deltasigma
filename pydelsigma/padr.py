import numpy as np

def padr(x, n, val=0.):
	"""y = padb(x, n, val)
	 Pad a matrix x on the right to length n with value val
	The empty matrix is assumed to be have 1 empty row
	"""
	if len(x.shape) == 1:
		xp = x.reshape((1, x.shape[0]))
	else:
		xp = x
	y = np.concatenate(
	                   (xp, 
	                    float(val)*np.ones((xp.shape[0], n - xp.shape[1]))
	                   ), axis=1
	                  )
	return y
	                  
def test_padr():
	tv = np.eye(15)
	tr = padr(tv, n=25, val=2)
	res = np.concatenate((tv, 2.*np.ones((15, 10))), axis=1)
	assert np.allclose(tr, res, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_padr()
