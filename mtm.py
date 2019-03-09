def ricker(f, dt):
    """
    Ricker wavelet of central frequency f
    
    IN    f : central frequency in Hz
          dt : time sampling interval in sec
          
    OUT   the Ricker wavelet
    
    Example
        wavelet = ricker(10, 0.002)
        
    Credit: Ke Chen, kechen@lbl.gov
    """
    
    tw = 2 * 1.35 * math.sqrt(6.0)/math.pi/f
    
    nw = 2 * math.floor(tw/dt/2) + 1

    nc = math.floor(nw/2)
    
    k = np.arange(1,nw+1)
    
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


from scipy import signal
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

# 1 - Read/synthetize data
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

tdelay = 0.0
s= shift_signal(d, tdelay)

plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Observed and synthetic data')
plt.show()

# 2 - Butterwoth bandpass filter data (optional)

# 3 - Windowing
n1 = 2800
n2 = 3200
d = d[n1:n2]
s = s[n1:n2]
t = np.arange(n1*dt,n2*dt,dt)
nt = len(t)
plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Windowed observed and synthetic data')
plt.show()

# 4 - compute_cc

# 5 - deconstruct_dat_cc (Apply CC -\delta T and -\delta A to the observed data prior to MTM)

# 6 - compute_average_error (sigma_t)

# 7 - FFT parameters
nextpow2 = math.ceil(math.log2(nt))
nfft = int(math.pow(2, nextpow2))
f0 = 0.0
df = 1.0/(nfft*dt)
dw = 2*math.pi*df
fnum = int(nfft/2) + 1
w = np.zeros(nfft)
w[0:fnum] = dw * np.arange(0,fnum)
w[fnum:nfft] = dw * np.arange(-fnum + 2, 0)

# 8 - DPSS
M = len(d)
NW = 2.5
Kmax = 5
ntaper = Kmax
dpss = signal.windows.dpss(M, NW, Kmax)
for i in range(0, Kmax):
    plt.plot(t,dpss[i,:])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('First five 2.5pi DPSS')
plt.show()

# 9 - Transfer function
wtr = 1.0e-10     #water level
T = np.zeros(fnum, dtype=complex)
A = np.zeros(fnum, dtype=complex)
B = np.zeros(fnum, dtype=complex)

for i in range(0, ntaper-1):
    #apply time domain taper
    dtp = d*dpss[i,:]
    stp = s*dpss[i,:]
    #apply FFT
    DTP = np.fft.fft(dtp,nfft)
    STP = np.fft.fft(stp,nfft)
    #
    for k in range(0,fnum):
        A[k] = A[k] + DTP[k] * np.conjugate(STP[i])
        B[k] = B[k] + STP[k] * np.conjugate(STP[i])

    plt.plot(t,dtp)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('DPSS tapered observed data')
plt.show()

for k in range(0, fnum):
    T[k] = A[k] / B[k]

#calculate phase
tau = np.zeros(fnum)

for k in range(0,fnum):
    tau[k] = math.atan2(T[k].imag, T[k].real)

faxis = np.arange(0,fnum)*df
plt.plot(faxis,tau)
plt.xlabel('Frequency')
plt.ylabel('Time delay')
plt.title('Frequency-dependent time delay')
plt.show()


