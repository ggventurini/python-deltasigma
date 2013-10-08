import numpy as np

def evalRPoly(roots, x, k=1):
	""" y = evalRPoly(roots, x, k=1)
	Compute the value of a polynomial which is given in terms of its roots.
	"""
	y = k * np.ones(x.shape)
	roots = roots[~np.isinf(roots)]        # remove roots at infinity
	for r in roots:
	    y = y*(x - r)
	return y
	
def test_evalRPoly():
	x = np.arange(1001)-500
	a = [1, 0, 1, 2]
	r1 = np.polyval(a, x)
	rts = np.roots(a)
	r2 = evalRPoly(rts, x)
	assert np.allclose(r1, r2, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_evalRPoly()
	
