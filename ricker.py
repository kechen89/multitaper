import math
import numpy as np

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
    
    k = np.arange(1,nw + 1)
    
    alpha = (k - nc - 1)*dt*f*math.pi
    
    beta = alpha*alpha
    
    w = (1.0 - 2.0*beta)*np.exp(-beta)
    
    return w
