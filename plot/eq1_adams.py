#!/usr/bin/env python

'''
    comments..
    
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#-----------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------
#a= np.fromfile('part.out',dtype=float ,count=-1 ,sep='')
#t_bwdE = np.loadtxt('bwdEuler001.out',usecols=(0,) )
#u_bwdE = np.loadtxt('bwdEuler001.out',usecols=(1,) )

#t_fwdE = np.loadtxt('fwdEuler001.out',usecols=(0,) )
#u_fwdE = np.loadtxt('fwdEuler001.out',usecols=(1,) )

real_t= np.loadtxt('realSol_001.out',usecols=(0,) )
real_y = np.loadtxt('realSol_001.out',usecols=(1,) )
 
ab4_t   = np.loadtxt('AB4_001.out',usecols=(0,) )
ab4_u   = np.loadtxt('AB4_001.out',usecols=(1,) )

 
am5_t   = np.loadtxt('AM5_001.out',usecols=(0,) )
am5_u   = np.loadtxt('AM5_001.out',usecols=(1,) )

ab5_t   = np.loadtxt('AB5_001.out',usecols=(0,) )
ab5_u   = np.loadtxt('AB5_001.out',usecols=(1,) )

 
#am_t   = np.loadtxt('AM3_001.out',usecols=(0,) )
#am3_u   = np.loadtxt('AM3_001.out',usecols=(1,) )

 
am4_t   = np.loadtxt('AM4_001.out',usecols=(0,) )
am4_u   = np.loadtxt('AM4_001.out',usecols=(1,) )

rk4_t   = np.loadtxt('RK4_001.out',usecols=(0,) )
rk4_u   = np.loadtxt('RK4_001.out',usecols=(1,) )









#b=

#i=0
#for row in a:
#    print(a[i] , '\t' , b[i])
#    i+=1


x = np.arange(0,2,0.25)
#y1 = fun_x(x)
#y2 = fun_x2(x)

#for i in range(len(x)):
#    print( x[i] , y1[i] )
majorLocator   = MultipleLocator(20)
majorFormatter = FormatStrFormatter('%f')
minorXLocator   = MultipleLocator(0.05)
minorYLocator   = MultipleLocator(0.025)






plt.rc('text', usetex=True )
#plt.rc('text', latex.preamble=r'\usepackage{babel}')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font',family='' ,size=20 )

fig,ax = plt.subplots(1,figsize=(12,6))
#plt.plot(real_t,real_y,'--k',linewidth=2,label=r'Analitical Sol.')
#plt.plot(t_bwdE,u_bwdE,'-.k',linewidth=2,label=r'Implicit Euler')
#plt.plot(t_fwdE,u_fwdE,'-.k',linewidth=2,label=r'Explicit Euler')

ax.set_xlim(0,2)


plt.plot(ab4_t,ab4_u,'gray',linewidth=2,label=r'Adams Bashforth 4 step')
plt.plot(ab5_t,ab5_u,'green',linewidth=2,label=r'Adams Bashforth 5 step')
plt.plot(am4_t,am4_u, 'c',linewidth=2,label=r'Adams Moulton 3 step' )
plt.plot(am5_t,am5_u,'r',linewidth=2,label=r'Adams Moulton 4 step')
plt.plot(rk4_t,rk4_u,'violet',linewidth=2,label=r'Runge Kutta 4 order')
plt.plot(real_t,real_y,'k-.',linewidth=2,label=r'Analitical Sol.')


plt.legend()
plt.title (r'$\dot{y} = 10(y-1) u$',fontsize=22)
plt.xlabel(r'$x$',fontsize=22)
plt.ylabel(r'$y$',fontsize=22)
#plt.xticks(range(0,2,0.5))
ax.xaxis.set_minor_locator(minorXLocator)
ax.yaxis.set_minor_locator(minorYLocator)

#plt.subplots_adjust(bottom=0.2)
#plt.margins()

#plt.grid(linestyle='dotted')
plt.grid(linestyle='--')
plt.tight_layout(0.5)
plt.savefig('case1.pdf')

plt.show()
'''
plt.figure(2, figsize=(10,6))
plt.plot(a,b,'r+')
plt.show()
'''
