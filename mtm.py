def ricker(f, dt):
    """
    Ricker wavelet
    """
    
    tw = 2 * 1.35 * math.sqrt(6.0)/math.pi/f
    
    nw = 2 * math.floor(tw/dt/2) + 1

    nc = math.floor(nw/2)
    
    k = np.arange(1,nw)
    
    alpha = (k - nc - 1)*dt*f*math.pi
    
    beta = alpha*alpha
    
    w = (1.0 - 2.0*beta)*np.exp(-beta)
    
    return w


import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

tmax = 12
dt = 0.002
nt = math.floor(tmax/dt) + 1
print(nt)
t = np.linspace(0,tmax,nt)
f = 1.0
wavelet = ricker(f, dt)

d = np.zeros(nt)
d[round(nt/2)] = 1.0
d = np.convolve(d,wavelet,'same')
#plt.plot(d)
#plt.show()

#time shift
nfft = int(math.pow(2, math.ceil(math.log2(nt))))
D = np.fft.fft(d,nfft)

S = np.zeros(nfft, dtype=complex)
delay = 1.0
for i in range(0, round(nfft/2)):
    w = 2.0*math.pi*i/nfft/dt
    Shift = cmath.exp(-1j*w*delay)
    S[i] = D[i]*Shift

# frequency domain symmetries
for i in range(1,round(nfft/2)-1):
    S[nfft - i] = np.conjugate(S[i])

s = np.fft.ifft(S)
s = np.real(s[0:nt])
print(len(s))
plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.show()



    


