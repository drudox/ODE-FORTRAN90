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
t_bwdE = np.loadtxt('bwdEuler001.out',usecols=(0,) )
u_bwdE = np.loadtxt('bwdEuler001.out',usecols=(1,) )

t_fwdE = np.loadtxt('fwdEuler001.out',usecols=(0,) )
u_fwdE = np.loadtxt('fwdEuler001.out',usecols=(1,) )

real_t= np.loadtxt('realSol_001.out',usecols=(0,) )
real_y = np.loadtxt('realSol_001.out',usecols=(1,) )
 
ab2_t   = np.loadtxt('AB2_001.out',usecols=(0,) )
ab2_u   = np.loadtxt('AB2_001.out',usecols=(1,) )

 
am2_t   = np.loadtxt('AM2_001.out',usecols=(0,) )
am2_u   = np.loadtxt('AM2_001.out',usecols=(1,) )





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


plt.plot(real_t,real_y,'k',linewidth=2,label=r'Analitical Sol.')
plt.plot(t_bwdE,u_bwdE,'r',linewidth=2,label=r'Implicit Euler')
plt.plot(t_fwdE,u_fwdE,'b',linewidth=2,label=r'Explicit Euler')
plt.plot(ab2_t,ab2_u,'g',linewidth=2,label=r'Adams Bashforth 2 step')



plt.legend()
plt.title (r'$\dot{y} = 10(y-1) u$',fontsize=22)
plt.xlabel(r'$x$',fontsize=22)
plt.ylabel(r'$y$',fontsize=22)
#plt.xticks(range(0,2,0.5))
ax.xaxis.set_minor_locator(minorXLocator)
ax.yaxis.set_minor_locator(minorYLocator)

#plt.subplots_adjust(bottom=0.2)
#plt.margins()

plt.grid()
plt.tight_layout(0.5)
plt.savefig('case1.pdf')

plt.show()
'''
plt.figure(2, figsize=(10,6))
plt.plot(a,b,'r+')
plt.show()
'''
