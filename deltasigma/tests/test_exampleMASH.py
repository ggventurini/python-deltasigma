from __future__ import division
import numpy as np
import deltasigma as ds
from scipy.signal import lti, ss2zpk, lfilter

def zpk_multiply(a, b):
    za, pa, ka = ds._utils._get_zpk(a)
    zb, pb, kb = ds._utils._get_zpk(b)
    pa = pa.tolist() if hasattr(pa, 'tolist') else pa
    pb = pb.tolist() if hasattr(pb, 'tolist') else pb
    za = za.tolist() if hasattr(za, 'tolist') else za
    zb = zb.tolist() if hasattr(zb, 'tolist') else zb
    return ds.cancelPZ((za+zb, pa+pb, ka*kb))

class testMultipleQ:
    def __init__(self):
        ABCD = [[1, 0, 0, 0, 1, -1, 0],
                [1, 1, 0, 0, 0, -2, 0],
                [0, 1, 1, 0, 0, 0, -1],
                [0, 0, 1, 1, 0, 0, -2],
                [0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0]]
        self.ABCD = np.array(ABCD, dtype=np.float_)
        self.nlev = [9, 9]
        k = [1., 1.]
        self.ntfs, self.stfs = ds.calculateTF(self.ABCD, k)

    def testMultipleQ1(self):
        """Test function for DS simulation with nq>1 1/2"""
        v1n = zpk_multiply(self.ntfs[0, 0], ([2, -1], [1, 0, 0, 0, 0]))
        v2n = zpk_multiply(self.ntfs[1, 0], ([1, 1], [0, 0], 1))
        # compute v1n/v2n and check that it is equal to -1
        res = zpk_multiply(v1n, (ds._utils._get_zpk(v2n)[1], ds._utils._get_zpk(v2n)[0], 1./ds._utils._get_zpk(v2n)[2]))
        assert int(ds.pretty_lti(res)) == -1

    def testMultipleQ2(self):
        """Test function for DS simulation with nq>1 2/2"""
        # filtering and simulation
        filtM1 = [0., 0., 0., 2., -1.]
        filtM2 = [1., -2., 1.]
        ntf_eq = zpk_multiply(self.ntfs[1, 1], self.ntfs[1, 1])
        M = self.nlev[0] - 1
        osr = 64
        f0 = 0.
        f1, f2 = ds.ds_f1f2(OSR=64, f0=0., complex_flag=False)
        delta = 2
        Amp = ds.undbv(-3) # Test tone amplitude, relative to full-scale.
        f = 0.3 # will be adjusted to a bin
        N = 2**12
        f1_bin = np.round(f1*N)
        f2_bin = np.round(f2*N)
        fin = np.round(((1 - f)/2*f1 + (f + 1)/2*f2) * N)
        # input sine
        t = np.arange(0, N).reshape((1, -1))
        u = Amp*M*np.cos((2*np.pi/N)*fin*t)
        vx, _, xmax, y = ds.simulateDSM(u, self.ABCD, nlev=self.nlev)
        # separate output #1 and output #2
        v1 = vx[0, :]
        v2 = vx[1, :]
        # filter and combine
        vf = lfilter(filtM1, [1.], v1) + lfilter(filtM2, [1.], v2)
        # compute the spectra
        window = ds.ds_hann(N)
        NBW = 1.5/N
        spec0 = np.fft.fft(vf*window)/(M*N/2)/ds.undbv(-6)
        spec1 = np.fft.fft(v1*window)/(M*N/2)/ds.undbv(-6)
        spec2 = np.fft.fft(v1*window)/(M*N/2)/ds.undbv(-6)
        freq = np.linspace(0, 0.5, N/2 + 1)

        # smooth, calculate the theorethical response and the SNR for VF
        spec0_smoothed = ds.circ_smooth(np.abs(spec0)**2., 16)
        Snn0 = np.abs(ds.evalTF(ntf_eq, np.exp(2j*np.pi*freq)))**2 * 2/12*(delta/M)**2
        snr0 = ds.calculateSNR(spec0[f1_bin:f2_bin + 1], fin - f1_bin)

        # smooth, calculate the theorethical response and the SNR for V1
        spec1_smoothed = ds.circ_smooth(np.abs(spec1)**2., 16)
        Snn1 = np.abs(ds.evalTF(self.ntfs[0, 0], np.exp(2j*np.pi*freq)))**2 * 2/12*(delta/M)**2
        snr1 = ds.calculateSNR(spec1[f1_bin:f2_bin + 1], fin - f1_bin)

        assert snr0 > 40
        assert snr1 > 40
        assert snr0-snr1 > 40

