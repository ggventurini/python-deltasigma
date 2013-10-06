import numpy as np

def zinc(f, m=64, n=1):
	""" mag = zinc(f, m=64, n=1)	
	Calculate the magnitude response of a cascade of n mth-order 
	comb filters at frequencies f.
	"""
	return np.fabs(np.sinc(m*f)/np.sinc(f) )**n
