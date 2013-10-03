import numpy as np

def thermometer(x, m):
	""" t = thermometer(x,m)	
	t is an m by length(x) matrix wherein the first
	x(i) components of column i are one
	"""
	t = np.zeros((m, len(x)))
	for i in range(len(x)):
	    t[:x[i], i] = np.ones((x[i], ))
	return t
	
def test():
	tv = np.arange(50)
	rm = np.zeros((70, 50))
	for i in range(50):
		rm[:i, i] = np.ones(rm[:i, i].shape)
	assert np.allclose(thermometer(tv, 70), rm, rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()

