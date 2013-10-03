import numpy as np

def undbv(x):
	""" y = undbv(x)	
	Convert x from dB to a voltage
	"""
	return 10.**(x/20.)

	
def test():
	assert np.allclose([undbv(53.05)], [449.26232467], rtol=1e-05, atol=1e-08)
	assert np.allclose([undbv(3)], [1.41253754462], rtol=1e-05, atol=1e-08)

if __name__ == '__main__':
	test()
