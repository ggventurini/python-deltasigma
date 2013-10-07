import numpy as np

def padl(x, n, val=0.):
	"""y = padb(x, n, val)
	 Pad a matrix x on the right to length n with value val
	The empty matrix is assumed to be have 1 empty row
	"""
	if len(x.shape) == 1:
		xp = x.reshape((1, x.shape[0]))
	else:
		xp = x
	y = np.concatenate(
	                   ( 
	                    float(val)*np.ones((max(1, x.shape[1]), n - x.shape[0])),
	                    xp
	                   ), axis=1
	                  )
	return y
	                  
def test_padl():
	tv = np.eye(15)
	tr = padl(tv, n=25, val=2)
	res = np.concatenate((2.*np.ones((15, 10)), tv), axis=1)
	assert np.allclose(tr, res, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_padl()
