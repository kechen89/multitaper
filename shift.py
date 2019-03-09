import math
import cmath
import numpy as np

def shift(d, dt, delay):
    """
        Shift a signal by time delay
    """
    nt = len(d)
    nextpow2 = math.ceil(math.log2(nt))
    nfft = int(math.pow(2, nextpow2))
    
    D = np.fft.fft(d,nfft)
    S = np.zeros(nfft, dtype=complex)
    
    for i in range(0, round(nfft/2) + 1):
        w = 2.0*math.pi*i/nfft/dt
        Shift = cmath.exp(-1j*w*delay)
        S[i] = D[i]*Shift
    
    # frequency domain symmetries
    for i in range(1,round(nfft/2)):
        S[nfft - i] = np.conjugate(S[i])

    s = np.fft.ifft(S)
    s = np.real(s[0:nt])
    return s
