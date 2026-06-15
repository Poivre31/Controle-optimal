import numpy as np
from matplotlib import pyplot as plt

t = np.linspace(0.5,3.5,1000)

dU = 1

Qopt = 2*t*dU
Qcste = dU*(2/(1-np.exp(-t))-1)
Qlin = dU*(2*((1+t)*np.exp(-t)+t*t/2-1)/(t-1+np.exp(-t))**2-1)

plt.plot(t,Qopt/Qopt)
plt.plot(t,Qcste/Qopt)
plt.plot(t,Qlin/Qopt)
plt.show()