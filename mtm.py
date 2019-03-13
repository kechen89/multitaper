from scipy import signal
import numpy as np
import math
import matplotlib.pyplot as plt
from ricker import ricker
from shift import shift

# 1 - Read/synthetize data
tmax = 12
dt = 0.002
nt = math.floor(tmax/dt) + 1
t = np.linspace(0,tmax,nt)
f = 1.0

wavelet = ricker(f, dt)
d = np.zeros(nt)
d[round(nt/2)] = 1.0
d = np.convolve(d,wavelet,'same')  #observed data

td = 0.002                      #time delay (if td > 0, syn arrive after obs)
s = shift(d, dt, td)            #synthetic data

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
d = d[n1:n2]    # d[n1] -> d[n2-1]
s = s[n1:n2]    # s[n1] -> s[n2-1]
t = np.arange(n1,n2)*dt # n1*dt -> (n2-1)*dt
nt = len(t)
plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Windowed observed and synthetic data')
plt.show()

# 4 - pre-process/time domain taper applied to windowed data
alpha = 10
nlen = len(s)
it_axis = np.arange(0,nlen)
cos_taper = 1.0 - np.cos(math.pi * it_axis / (nlen - 1)) ** alpha
d = d * cos_taper
s = s * cos_taper
plt.plot(t, cos_taper)
plt.title('Cosine taper')
plt.show()

plt.plot(t,d,'r')
plt.plot(t,s,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Tapered observed and synthetic data')
plt.show()

# 5 - compute_cc
TSHIFT_MIN = -5.0
TSHIFT_MAX = 5.0
cc = np.correlate(d, s, "same")
#norm_s = math.sqrt(sum(s*s))
#norm_d = math.sqrt(sum(d*d))
#norm = norm_s * norm_d
#cc = cc/norm
ishift = np.argmax(cc) - int(nlen/2)
tshift = ishift*dt
dlnA = 0.5 * math.log(sum(d*d) / sum(s*s) )
print('Time shift measured by cross-correlation:', tshift, 's')
print('Amplitude difference measured by cross-correlation:', dlnA)

# 6 - deconstruct_dat_cc (Apply CC -\delta T and -\delta A to the observed data prior to MTM)
d_dec = np.zeros(nlen)
for i in range(0, nlen):
    if (i + ishift) >= 0 and (i + ishift) <= nlen - 1:
        d_dec[i] = d[i + ishift]
print(d_dec)
#if (ishift < 0) d_dec(1:-ishift+1) = dat_dtw_cc(-ishift+2)
#if (ishift > 0) dat_dtw_cc(nlen-ishift:nlen) = dat_dtw_cc(nlen-ishift-1)

plt.plot(t,s,'r')
plt.plot(t,d_dec,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Deconstruct data')
plt.show()

# 7 - compute_average_error (sigma_t)
#function
# compute_average_error

# 8 - FFT parameters
nextpow2 = math.ceil(math.log2(nt))
nfft = int(math.pow(2, nextpow2))
f0 = 0.0
df = 1.0/(nfft*dt)
dw = 2*math.pi*df
fnum = int(nfft/2) + 1

w = np.zeros(nfft)    # angular frequency vector
w[0:fnum] = dw * np.arange(0,fnum)
w[fnum:nfft] = dw * np.arange(-fnum + 2, 0)

# 9 - DPSS
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

# 10 - Transfer function
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

for k in range(1,fnum):
    tau[k] = math.atan2(T[k].imag, T[k].real)

faxis = np.arange(0,fnum)*df
plt.plot(faxis,tau)
plt.xlabel('Frequency')
plt.ylabel('Time delay')
plt.title('Frequency-dependent time delay')
plt.show()


