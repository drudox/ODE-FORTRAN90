
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
from matplotlib import colors as mcolors
#import 
plt.style.use(['marco.mplstyle'])
#plt.style.use(['marco.mplstyle'])

x = np.genfromtxt('datiFWE_1.dat',usecols=0,)
x2 = np.genfromtxt('datiBWE_1.dat',usecols=0,)
x3 = np.genfromtxt('realdata1.dat',usecols=0,)
x4 = np.genfromtxt('datiRK4_1.dat',usecols=0,)
y = np.genfromtxt('datiFWE_1.dat',usecols=1,)
y2 = np.genfromtxt('datiBWE_1.dat',usecols=1,)
y3 = np.genfromtxt('realdata1.dat',usecols=1,)
y4 = np.genfromtxt('datiRK4_1.dat',usecols=1,)

fig, ax = plt.subplots(1, figsize=(6, 3))


ax.plot(x, y, color='navy' ,linewidth=2, label=r'FWD Euler')
ax.plot(x2, y2, color='red' ,linewidth=2, label=r'BWD Euler')
ax.plot(x3, y3, color='k' ,linewidth=1, label=r'Analitical')
ax.plot(x4, y4, color='g', linestyle='-.' ,linewidth=2, label=r'RK-4th')
#ax.plot(x, y2, 'k'  , linewidth=2, label=r'$e^{-x^2}$')
#ax.plot(x, y1, color='maroon' ,linewidth=2, label=r'$e^{-x}$')
ax.locator_params(axis='y', tight=True, nbins=5)
ax.locator_params(axis='x',  tight=True , nbins=20)

ax.minorticks_on()
ax.legend()
plt.tight_layout(0.5)
#plt.show()
plt.title(r'yp=-10(x-1)t') 
plt.savefig('euler1.pdf')

