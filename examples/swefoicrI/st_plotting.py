# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 05:31:39 2014

@author: willem
"""

#import pylab as pl
#import numpy as np
import matplotlib.pyplot as plt

f = open('../sinybed/out.160')
a = f.read().split()
vardims = int(a[0])
something = int(a[1])
griddim = int(a[2])
x = [float(d) for d in a[3:3+griddim*2:2]]
n = 1
y1 = [float(d) for d in a[4+griddim*(2*n-2):4+griddim*2*n:2]]
n = 2
y2 = [float(d) for d in a[4+griddim*(2*n-2):4+griddim*2*n:2]]
n = 3
y3 = [float(d) for d in a[4+griddim*(2*n-2):4+griddim*2*n:2]]
plt.plot(x,y3,x,y2,x,[y1[k] + y3[k] for k in range(len(y1))])
plt.savefig('../sinybed/test.png')
plt.show()



#x = np.linspace(0, 10)
#line, = plt.plot(x, np.sin(x), '--', linewidth=2)
#
#dashes = [10, 5, 100, 5] # 10 points on, 5 off, 100 on, 5 off
#line.set_dashes(dashes)
#
#
#plt.show()