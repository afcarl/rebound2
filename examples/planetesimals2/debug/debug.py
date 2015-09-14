#purpose of this macro is to compare different Hybrid runs from each other, and plot the x,y,z, differences.

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names=['time (years)','x','y','z','vx','vy','vz']
colors=['b','g','m','r','c','y']
bodynames=['star','planet 1','planet 2','particle']

#file_name1 = 'HYB_sept14_single.txt'
file_name1 = 'HYB_sept14_par11.txt'


N_active = 4

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=0
arg2=int(sys.argv[1])
arg3 = -1
arg4 = 0
if len(sys.argv) == 3:
    arg3 = int(sys.argv[2])
if len(sys.argv) == 4:
    arg4 = int(sys.argv[3])

fos = open(file_name1, 'r')
data = np.loadtxt(fos, delimiter=',')

for i in xrange(0,1):
    p=data[i::N_active]
    plt.plot(p[arg4:arg3,arg1], p[arg4:arg3,arg2], 'o'+colors[i], marker='o', markersize=5-i,markeredgecolor='none',label=bodynames[i] )

#plt.xlim([33,45])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.title(file_name1)
plt.show()
