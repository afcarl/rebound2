#purpose of this macro is to compare different Hybrid runs from each other, and plot the x,y,z, differences.

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names=['time (years)','dx','dy','dz','dvx','dvy','dvz']
colors=['b','g','m','r','c','y']
bodynames=['star','planet 1','planet 2','particle']

#file_name1 = 'HYB_sept14_single.txt'

#file_name1 = 'HYB_sept15_par13_single.txt'
#file_name2 = 'HYB_sept15_par13_group.txt'

#file_name1 = 'IAS_sept16_single.txt'
#file_name2 = 'IAS_sept16.txt'
#N_active = 4

file_name1 = 'IAS_sept16_pforces_single.txt'
file_name2 = 'IAS_sept16_pforces.txt'
N_active = 3

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=0
arg2=int(sys.argv[1])
arg3 = -1
arg4 = 0
if len(sys.argv) >= 3:
    arg3 = int(sys.argv[2])
if len(sys.argv) >= 4:
    arg4 = int(sys.argv[3])

fos = open(file_name1, 'r')
data = np.loadtxt(fos, delimiter=',')
fos2 = open(file_name2, 'r')
data2 = np.loadtxt(fos2, delimiter=',')

for i in xrange(0,N_active):
    p=data[i::N_active]
    q=data2[i::N_active]
    plt.plot(p[arg4:arg3,arg1], p[arg4:arg3,arg2] - q[arg4:arg3,arg2], 'o'+colors[i], marker='o', markersize=2,markeredgecolor='none',label=bodynames[i] )
    #plt.plot(p[arg4:arg3,arg1], q[arg4:arg3,arg2], color=colors[i], linestyle='--',label=bodynames[i] )

#plt.xlim([33,45])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.title(file_name1)
plt.show()
