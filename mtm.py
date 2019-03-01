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
import matplotlib.pyplot as plt

tmax = 12
dt = 0.002
nt = math.floor(tmax/dt) + 1
t = np.linspace(0,tmax,nt)
f = 1.0
w = ricker(f, dt)

d = np.zeros(nt)
d[round(nt/2)] = 1.0
d = np.convolve(d,w,'same')
plt.plot(d)
plt.show()

#time shift





    


