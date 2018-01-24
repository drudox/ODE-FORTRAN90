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
t_rkf = np.loadtxt('../result/rungekutta/RKFehlberg_001.out',usecols=(0,) )
u_rkf = np.loadtxt('../result/rungekutta/RKFehlberg_001.out',usecols=(1,) )

t_rk4 = np.loadtxt('../result/rungekutta/RK4_001.out',usecols=(0,) )
u_rk4 = np.loadtxt('../result/rungekutta/RK4_001.out',usecols=(1,) )



t_ab2 = np.loadtxt('../result/adams/AB2_001.out',usecols=(0,) )
u_ab2 = np.loadtxt('../result/adams/AB2_001.out',usecols=(1,) )


t_am2 = np.loadtxt('../result/adams/AM2_001.out',usecols=(0,) )
u_am2 = np.loadtxt('../result/adams/AM2_001.out',usecols=(1,) )


t_me = np.loadtxt('../result/rungekutta/Heun_001.out',usecols=(0,) )
u_me = np.loadtxt('../result/rungekutta/Heun_001.out',usecols=(1,) )

t_lf = np.loadtxt('../result/multistep/LeapFrog_001.out',usecols=(0,) )
u_lf = np.loadtxt('../result/multistep/LeapFrog_001.out',usecols=(1,) )


real_t= np.loadtxt('../realSol_001.out',usecols=(0,) )
real_y = np.loadtxt('../realSol_001.out',usecols=(1,) )
 

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

fig,ax = plt.subplots(1,figsize=(9,6))
#plt.plot(real_t,real_y,'--k',linewidth=2,label=r'Analitical Sol.')
#plt.plot(t_bwdE,u_bwdE,'-.k',linewidth=2,label=r'Implicit Euler')
#plt.plot(t_fwdE,u_fwdE,'-.k',linewidth=2,label=r'Explicit Euler')

ax.set_xlim(0,2)


#plt.plot(real_t,real_y,'k',linewidth=1,label=r'Analitical Sol.')
#plt.plot(t_cn,u_cn,'k',linewidth=1,linestyle=':',label=r'Crank Nicolson')
#plt.plot(t_me,u_me,'k',linewidth=1,linestyle='-.',label=r'Modified Euler')
#plt.plot(t_lf,u_lf,'k',linewidth=1,linestyle='--',label=r'Leap Frog')

plt.plot(real_t,real_y,'k',linewidth=4,label=r'Analitical Sol.')
plt.plot(t_rkf,u_rkf,'b',linewidth=2,label=r'Runge Kutta Fehlberg 4-5th')
plt.plot(t_rkf,u_rkf,'b+',linewidth=2,markersize=14) #,label=r'e Kutta Fehlberg 4-5th')
plt.plot(t_rk4,u_rk4,'r',linewidth=1,linestyle='-',label=r'Runge Kutta 4th')

#plt.plot(t_lf,u_lf,'g',linewidth=2,linestyle='-',label=r'Leap Frog')



plt.legend()
plt.title (r'$\dot{y} = 10(t-1) y$',fontsize=22)
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
plt.savefig('Eq1_RungeKutta.pdf')

plt.show()
'''
plt.show()
'''
