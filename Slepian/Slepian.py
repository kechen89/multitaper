from scipy import signal
import numpy as np
import math
import matplotlib.pyplot as plt
import struct

# 1 - Read/synthetize data

s = np.fromfile('DPSS.bin', dtype=float)
s = np.reshape(s,(4,64))

for i in range(0, 4):
    plt.plot(s[i,:])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('DPSS')
plt.show()


# 9 - DPSS
dpss = signal.windows.dpss(64, 4, 4)

for i in range(0, 4):
    plt.plot(dpss[i,:])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('DPSS')
plt.show()


diff = s-dpss
for i in range(0, 4):
    plt.plot(diff[i,:])
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Difference')
plt.show()
