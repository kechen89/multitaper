from scipy import signal
import numpy as np
import math
import matplotlib.pyplot as plt
from utils import ricker
from utils import shift
from obspy import read

# 1 - Read/synthetize data

# Synthetize data
#tmax = 12
#dt = 0.002
#nt = math.floor(tmax/dt) + 1
#t = np.arange(0,nt)*dt
#f = 1.0

#wavelet = ricker(f, dt)
#d = np.zeros(nt)
#d[round(nt/2)] = 1.0
#d = np.convolve(d,wavelet,'same')  #observed data

#td = 0.1                        #time delay (if td > 0, syn arrive after obs)
#s = shift(d, dt, td)            #synthetic data

# Read data
st = read('trival_data.sac',debug_headers=True)
d = st[0].data
st = read('trival_syn.sac',debug_headers=True)
s = st[0].data
dt = st[0].stats.delta
nt = st[0].stats.npts
t = np.arange(0,nt)*dt

plt.plot(t,d,'b')
plt.plot(t,s,'r')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Observed (blue) and synthetic data (red)')
plt.show()
plt.magnitude_spectrum(d,Fs=1/dt)
plt.title('Frequency Spectrum of data')
plt.show()

# 2 - Butterwoth bandpass filter data (optional)

# 3 - Windowing
n1 = 0
n2 = nt

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
cos_taper = 1.0 - (np.cos(math.pi * it_axis / (nlen - 1))) ** alpha
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

# 5 - Compute_cc
cc = np.correlate(d, s, "same")       # cc is not normalized
ishift = np.argmax(cc) - int(nlen/2)
tshift = ishift*dt
dlnA = 0.5 * math.log(sum(d*d) / sum(s*s) )
print('Time shift measured by cc:', tshift, 's')
print('Amplitude difference measured by cc:', dlnA)

# 6 - Deconstruct_dat_cc (Apply CC -\delta T and -\delta A to the observed data prior to MTM)
# Why?
d_dec = np.zeros(nlen)
for i in range(0, nlen):
    if (i + ishift) > 0 and (i + ishift) < nlen - 1:
        d_dec[i] = d[i + ishift]
if ishift < 0:
    d_dec[0:-ishift] = d_dec[-ishift+1]
if ishift > 0:
    d_dec[nlen-ishift-1:nlen-1] = d_dec[nlen-ishift-2]

d_dec = d_dec * math.exp(-dlnA)

plt.plot(t,s,'r')
plt.plot(t,d_dec,'b')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Synthetic data and Deconstruct data')
plt.show()

# 7 - Compute the estimated uncertainty for the cross-correlation measurement
# based on integrated waveform difference between the data and reconstructed synthetics
# compute_average_error

# 8 - FFT parameters
nextpow2 = math.ceil(math.log2(nt))
nfft = 2 * int(math.pow(2, nextpow2))
f0 = 0.0
df = 1.0/(nfft * dt)
dw = 2 * math.pi * df
fnum = int(nfft/2) + 1

fvec = np.zeros(nfft)    # angular frequency vector
fvec[0:fnum] = df * np.arange(0,fnum)    # positive frequency
fvec[fnum:nfft] = df * np.arange(-fnum + 2, 0)  # negative frequency

wvec = np.zeros(nfft)    # angular frequency vector
wvec[0:fnum] = dw * np.arange(0,fnum)    # positive frequency
wvec[fnum:nfft] = dw * np.arange(-fnum + 2, 0)  # negative frequency

# Estimate stopping frequency
WTR = 0.02
ampmax = 0.0
k_amp_max = 0
S = np.fft.fft(s,nfft)
D = np.fft.fft(d,nfft)

for k in range(0, fnum):
    if  abs(S[k]) > ampmax:
        ampmax =  abs(S[k])
        k_amp_max = k

wtr_level = ampmax * WTR    # water level value to stablize

fmax = fnum
fmax_stop = 0

for k in range(0, fnum):
    if abs(S[k]) <= wtr_level and fmax_stop == 0 and k > k_amp_max:
        fmax_stop = 1
        fmax = k
    if abs(S[k]) >= 10.0 * wtr_level and fmax_stop == 1 and k > k_amp_max:
        fmax_stop = 0
        fmax = k

print('Stopping frequency index',fmax)
print('Stopping frequency',fmax*df, 'Hz')

#Spectrum of synthetic data
S_SPEC = np.abs(S)
D_SPEC = np.abs(D)
plt.plot(wvec[0:fmax],S_SPEC[0:fmax])
plt.plot(wvec[0:fmax],D_SPEC[0:fmax])
plt.title('Spectrum of synthetic and observed data')
axes = plt.gca()
axes.set_xlim([0,1.2])
plt.show()

# 9 - DPSS
M = len(d)
NW = 0.5
Kmax = int(2.0 * NW)
ntaper = Kmax
dpss = signal.windows.dpss(M, NW, Kmax)

for i in range(0, Kmax):
    plt.plot(t,dpss[i,:])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('DPSS')
plt.show()

# 10 - Transfer function
T = np.zeros(fnum, dtype=complex)
A = np.zeros(fnum, dtype=complex)
B = np.zeros(fnum, dtype=complex)

for i in range(0, ntaper):
    #apply time domain taper
    #dtp = d_dec * dpss[i,:]
    dtp = d * dpss[i,:]
    stp = s * dpss[i,:]
    #apply FFT
    DTP = np.fft.fft(dtp,nfft)
    STP = np.fft.fft(stp,nfft)
    #
    for k in range(0,fnum):      # sum over taper
        A[k] = A[k] + DTP[k] * np.conjugate(STP[k])
        B[k] = B[k] + STP[k] * np.conjugate(STP[k])

    plt.plot(t,dtp)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('DPSS tapered observed data')
plt.show()

# water level
wtr = 1e-2     #water level
ampmax = 0.0
i_amp_max = 1

for k in range(0, fnum):
    if  abs(B[k]) > ampmax:
        ampmax =  abs(B[k])
        i_amp_max = k
epsilon = ampmax * wtr**2    # water level value to stablize
print(epsilon)

for k in range(0, fnum):
    if abs(B[k]) > epsilon:
        T[k] = A[k] / B[k]
    else:
        T[k] = A[k] / (B[k] + epsilon)

# 11 - Calculate phase
phase = np.zeros(fmax)
dtau = np.zeros(fmax)

for k in range(0,fmax):
    phase[k] = math.atan2(T[k].imag, T[k].real)

waxis = np.arange(0,fmax)*dw
plt.plot(waxis,phase)
plt.title('Phase')
plt.xlabel('Angular frequency')
plt.show()

phase1 = np.unwrap(phase)
plt.plot(waxis,phase1)
plt.title('Unwrapped Phase')
plt.xlabel('Angular frequency')
plt.show()

# 12 - Phase correction and Calculate time delay
# phase correction parameters, between (PI, 2PI), use a higher value for conservative phase wrapping
PHASE_STEP = 1.5 * math.pi

#dtau[0] = tshift
for k in range(0, fmax):
    if k > 0 and k < fmax - 1 :
        
        smth = phase[k + 1] + phase[k - 1] - 2.0 * phase[k]
        smth1 = (phase[k + 1] + 2.0*math.pi) + phase[k - 1] - 2.0 * phase[k]
        smth2 = (phase[k + 1] - 2.0*math.pi) + phase[k - 1] - 2.0 * phase[k]
        
        if abs(smth1) < abs(smth) and abs(smth1) < abs(smth2) and abs(phase[k] - phase[k+1]) > PHASE_STEP:
            for j in range(k + 1, fmax):
                phase[j] = phase[j] + 2.0*math.pi
        if abs(smth2) < abs(smth) and abs(smth2) < abs(smth1) and abs(phase[k] - phase[k+1]) > PHASE_STEP:
            for j in range(k+1, fmax):
                phase[j] = phase[j] - 2.0*math.pi
    if k > 0:
        dtau[k] = (-1.0/wvec[k]) * phase[k] #+ tshift

waxis = np.arange(0,fmax)*dw
plt.plot(waxis,phase)
plt.xlabel('Angular Frequency')
plt.ylabel('Phase')
plt.title('Frequency-dependent Phase')
plt.show()

plt.plot(waxis,dtau)
axes = plt.gca()
axes.set_ylim([0,4])
plt.xlabel('Angular Frequency')
plt.ylabel('Time delay')
plt.title('Time delay')
plt.show()

# estimate error using CC

# 13 - Adjust frequency range for MT measurements
i_left = 0
i_right = fmax
fstart = max(0, 1.5/(nt*dt))
fend = 1.0/(2.0*dt)

i_left_old = i_left
i_right_old = i_right

for j in range(i_left_old, i_right_old):
    if fvec[j] > fstart :
        i_left = j
        break

for j in range(i_left_old, i_right_old):
    if fvec[j] > fend :
        i_right = j
        break

fstart = i_left * df
fend = i_right * df

print('New frequency band:', fstart, fend)
waxis = np.arange(i_left, i_right) * dw
dtau_new = dtau[i_left:i_right]
plt.plot(waxis,dtau_new)
axes = plt.gca()
axes.set_ylim([0,4])
plt.xlabel('Frequency')
plt.ylabel('Time delay')
plt.title('Time delay')
plt.show()

# 14 - Misfit function (integration over frequency)
w_taper = np.zeros(fnum)

for k in range(0, fmax):
    w_taper[k] = 1.0 - math.cos(math.pi * k/(fmax - 1))**alpha

misfit = 0.5 * 2.0 * df * np.sum((dtau[0:fmax])**2 * w_taper[0:fmax])

print('misfit:',misfit)





