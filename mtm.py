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

def shift_signal(d, delay):
    """
    Shift a signal by time delay
    """
    nextpow2 = math.ceil(math.log2(nt))
    nfft = int(math.pow(2, nextpow2))
    D = np.fft.fft(d,nfft)
    S = np.zeros(nfft, dtype=complex)

    for i in range(0, round(nfft/2)):
        w = 2.0*math.pi*i/nfft/dt
        Shift = cmath.exp(-1j*w*delay)
        S[i] = D[i]*Shift

    # frequency domain symmetries
    for i in range(1,round(nfft/2)-1):
        S[nfft - i] = np.conjugate(S[i])

    s = np.fft.ifft(S)
    s = np.real(s[0:nt])
    return s


from scipy import signal
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

s= shift_signal(d, 1.0)

plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Observed and synthetic data')
plt.show()

#window signal 4-9s
d = d[2000:4500]
s = s[2000:4500]
t = np.arange(2000*dt,4500*dt,dt)
print(len(t))
plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Windowed observed and synthetic data')
plt.show()

#DPSS
M = len(d)
NW = 2.5
Kmax = 5
dpss = signal.windows.dpss(M, NW, Kmax)
for i in range(0, Kmax):
    plt.plot(t,dpss[i,:])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('First five 2.5pi DPSS')
plt.show()

#Multitaper misfit


