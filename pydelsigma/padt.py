import numpy as np

def padt(x, n, val=0.):
	"""y = padb(x, n, val)
	Pad a matrix x on the top to length n with value val(0)
	The empty matrix is assumed to be have 1 empty column.
	"""
	if len(x.shape) == 1:
		xp = x.reshape((x.shape[0], 1))
	else:
		xp = x
	y = np.concatenate(
	                   (float(val)*np.ones((n - xp.shape[0], xp.shape[1])),
	                    xp
	                   ), axis=0
	                  )
	return y
	                  
def test_padt():
	tv = np.eye(15)
	tr = padt(tv, n=25, val=2)
	res = np.concatenate((2.*np.ones((10, 15)), tv), axis=0)
	assert np.allclose(tr, res, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_padt()
