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
t_bwdE = np.loadtxt('../result/euler/bwdEuler_002.out',usecols=(0,) )
u_bwdE = np.loadtxt('../result/euler/bwdEuler_002.out',usecols=(1,) )

t_fwdE = np.loadtxt('../result/euler/fwdEuler_002.out',usecols=(0,) )
u_fwdE = np.loadtxt('../result/euler/fwdEuler_002.out',usecols=(1,) )

real_t= np.loadtxt('../realSol_002.out',usecols=(0,) )
real_y = np.loadtxt('../realSol_002.out',usecols=(1,) )
 

#b=

#i=0
#for row in a:
#    print(a[i] , '\t' , b[i])
#    i+=1


x = np.arange(0,2.5)
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

fig,ax = plt.subplots(1,figsize=(9,6))
#plt.plot(real_t,real_y,'--k',linewidth=2,label=r'Analitical Sol.')
#plt.plot(t_bwdE,u_bwdE,'-.k',linewidth=2,label=r'Implicit Euler')
#plt.plot(t_fwdE,u_fwdE,'-.k',linewidth=2,label=r'Explicit Euler')

ax.set_xlim(0,2.5)


plt.plot(real_t,real_y,'k',linewidth=2,linestyle='-',label=r'Analitical Sol.')
plt.plot(t_bwdE,u_bwdE,'k',linewidth=2,linestyle=':',label=r'Implicit Euler')
plt.plot(t_fwdE,u_fwdE,'k',linewidth=2,linestyle='-.',label=r'Explicit Euler')



plt.legend()
plt.title (r'$\dot{y} = -20y+20\sin(t)+\cos(t) $',fontsize=22)
plt.xlabel(r'$x$',fontsize=22)
plt.ylabel(r'$y$',fontsize=22)
#plt.xticks(range(0,2,0.5))
ax.xaxis.set_minor_locator(minorXLocator)
ax.yaxis.set_minor_locator(minorYLocator)

#plt.subplots_adjust(bottom=0.2)
#plt.margins()

#plt.grid()
plt.grid(linestyle='dotted')
plt.tight_layout(0.5)
plt.savefig('Euler_eq2.pdf')

plt.show()
'''
plt.figure(2, figsize=(10,6))
plt.plot(a,b,'r+')
plt.show()
'''
