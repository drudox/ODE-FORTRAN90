#!/usr/bin/env python

'''
    comments..
    
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#-----------------------------------------------------------------------------------------------------

def fun_x(x):
    y= np.exp(-x)
    return y

def fun_x2(x):
    y= np.exp(-x**2)
    return y




#-----------------------------------------------------------------------------------------------------
#a= np.fromfile('part.out',dtype=float ,count=-1 ,sep='')
a= np.loadtxt('part.out',usecols=(0,) )
b= np.loadtxt('part.out',usecols=(1,) )

#b=

i=0
for row in a:
    print(a[i] , '\t' , b[i])
    i+=1


x = np.arange(0,10,0.1)
y1 = fun_x(x)
y2 = fun_x2(x)

#for i in range(len(x)):
#    print( x[i] , y1[i] )
majorLocator   = MultipleLocator(20)
majorFormatter = FormatStrFormatter('%d')
minorXLocator   = MultipleLocator(0.25)
minorYLocator   = MultipleLocator(0.025)






plt.rc('text', usetex=True )
#plt.rc('text', latex.preamble=r'\usepackage{babel}')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('font',family='' ,size=20 )

fig,ax = plt.subplots(1,figsize=(9,6))
plt.plot(x,y1,'--k',linewidth=2,label=r'$e^{-x}$')
plt.plot(x,y2,'-.k',linewidth=2,label=r'$e^{-x^2}$')

plt.legend()
plt.title (r'$\text{Order of} \iint \exp$',fontsize=22)
plt.xlabel(r'$x$',fontsize=22)
plt.ylabel(r'$y$',fontsize=22)
plt.xticks(range(0,10,1))
ax.xaxis.set_minor_locator(minorXLocator)
ax.yaxis.set_minor_locator(minorYLocator)

#plt.subplots_adjust(bottom=0.2)
#plt.margins()


plt.tight_layout(0.5)
plt.savefig('graphic.pdf')

plt.show()
'''
plt.figure(2, figsize=(10,6))
plt.plot(a,b,'r+')
plt.show()
'''
