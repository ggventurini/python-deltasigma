import numpy as np
import control
from control.xferfcn import TransferFunction
from control.statesp import StateSpace
from pydelsigma import evalRPoly

def evalTF(tf, z):
	"""h = evalTF(tf, z)
	Evaluates the rational function described by tf
	at the point(s) given in the z vector.
	
	TF must be either a TransferFunction or an LTI object or any object having the attributes:
	 zeros, poles, k	if form == 'zp'
	or:
	 num, den		if form == 'coeff'
	"""
	#In Matlab 5, the ss/freqresp() function does nearly the same thing.
	
	# in our case a transfer function is a 'TransferFunction' not a zpk
	if (hasattr(tf, 'inputs') and not tf.inputs == 1) or \
	   (hasattr(tf, 'outputs') and not tf.outputs == 1):
			raise TypeError, "Only SISO transfer functions can be evaluated."
	if hasattr(tf, 'num') and hasattr(tf, 'den'):
		filt = hasattr(tf, 'outputs')
		num = tf.num[0][0] if filt else tf.num
		den = tf.den[0][0] if filt else tf.den
		h = np.polyval(num, z) / np.polyval(den, z)
	elif (hasattr(tf, 'zeros') and hasattr(tf, 'poles')) or \
	   (hasattr(tf, 'zero') and hasattr(tf, 'pole')):
		# LTI objects have poles and zeros, 
		# TransferFunction-s have pole() and zero()
	   	zeros = tf.zeros if hasattr(tf, 'zeros') else tf.zero()
	   	poles = tf.poles if hasattr(tf, 'poles') else tf.pole()
		if hasattr(tf, 'k'):
			k = tf.k
		elif hasattr(tf, 'gain'):
			k = tf.gain  
		elif hasattr(tf, 'returnScipySignalLti'): 
			k = np.array(tf.returnScipySignalLti()[0][0].gain)
		h = k * evalRPoly(zeros, z) / evalRPoly(poles, z)
	elif hasattr(tf, 'form') and tf.form == 'zp':
		h = tf.k * evalRPoly(tf.zeros, z) / evalRPoly(tf.poles, z)
	elif hasattr(tf, 'form') and tf.form == 'coeff':
		h = np.polyval(tf.num, z) / np.polyval(tf.den, z);
	elif hasattr(tf, 'form'):
		raise ValueError, '%s: Unknown form: %s' % (__name__, tf.form)
	else:
		raise TypeError, '%s: Unknown transfer function %s' % (__name__, str(tf))
	return h
	
def test_evalTF():
	from control.matlab import tf, tf2zpk
	from pydelsigma.utils import empty
	num, den = np.poly([3, 0.3, 1]), np.poly([2, 0.5, .25])
	H = tf(num, den, 1)
	tstr1 = empty()
	tstr1.form, tstr1.num, tstr1.den = 'coeff', num, den
	tstr2 = empty()
	tstr2.form = 'zp'
	tstr2.zeros, tstr2.poles, tstr2.k = tf2zpk(num, den)
	z = np.exp(1j*np.linspace(0, 2*np.pi, num=129, endpoint=True))
	h1 = evalTF(tstr1, z)
	h2 = evalTF(tstr2, z)
	h3 = evalTF(H, z)
	assert np.allclose(np.abs(h1), np.abs(h2), atol=1e-8, rtol=1e-5)
	assert np.allclose(h1, h3, atol=1e-8, rtol=1e-5)

if __name__ == '__main__':
	test_evalTF()
	
