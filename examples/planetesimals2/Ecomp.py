#purpose of this macro is to compare the HYBRID->IAS integrator energies to the HYBRID integrator energies. A debugging mechanism

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names=['time (years)','Semi-Major Axis (AU)','DEccentricity','D(Ei - E0) / E0','D(Ki - K0) / K0','D(Ui - U0) / U0','DTotal Ang. Mom.','Dplanet-star distance']
colors=['b','g','m','r','c','y']

file_name1 = 'Ecomp/planet_HYBRID_5_0.txt'
file_name2 = 'Ecomp/planet_HYBRIDIAS_5_0.txt'

N_active = 2

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
fos2 = open(file_name2, 'r')
data2 = np.loadtxt(fos2, delimiter=',')

for i in xrange(0,N_active):
    p=data[i::N_active]
    q=data2[i::N_active]
    plt.plot(p[arg4:arg3,arg1], p[arg4:arg3,arg2]-q[arg4:arg3,arg2], 'o'+colors[i], marker='o', markersize=2,markeredgecolor='none',label='planet '+str(i)+' positive' )
    plt.plot(p[arg4:arg3,arg1], q[arg4:arg3,arg2]-p[arg4:arg3,arg2], 'o', color='red', marker='o', markersize=2,markeredgecolor='none',label='planet '+str(i)+' negative' )
    if arg2 == 3 or arg2==4 or arg2==5 or arg2==6:
        plt.yscale('log')
        break

plt.xlim([p[arg4,0],p[arg3,0]])
#plt.xlim([33,45])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.title('Difference of HYBRID - HYBRIDIAS')
plt.show()
