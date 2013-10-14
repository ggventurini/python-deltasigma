import numpy as np

def axisLabels(ran, incr):
	"""function s = axisLabels(ran, incr)
	range is an array containing the axis points (floats)
	incr may be:
	* an int, the function returns an array of strings corresponding to:
	each element of range[0]:range[-1]:incr formatted as '%g'
	* a list, the function returns an array of strings corresponding to:
	each element of incr[1]:range[-1]:incr[0] formatted as '%g'
	
	Note: all elements in ran less than 1e-6 are rounded down to 0.
	"""
	ran[np.abs(ran)<1e-6] = 0
	s = []
	if not hasattr(incr, '__len__'):
		incr = int(incr)
		first = 0
	elif len(incr) == 2:
		first = incr[1]
		incr = incr[0]
	else:
		raise ValueError, "Unrecognised incr: "+str(incr)
	for i in range(first, len(ran), incr):
		s += ['%g' % ran[i]]
	return s
	
def test_axisLabels():
	ran = np.arange(100)
	ss = axisLabels(ran, incr=10)
	r = ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90']
	assert r == ss
	
if __name__ == '__main__':
	test_axisLabels()
