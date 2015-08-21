import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names=['time (years)','Semi-Major Axis (AU)','Eccentricity','(Ei - E0) / E0','Total Ang. Mom.']
colors=['b','g','m','r','c','y']
N=1

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=int(sys.argv[2])
arg2=int(sys.argv[3])
arg3 = -1
arg4 = 0
if len(sys.argv) == 5:
    arg3 = int(sys.argv[4])
elif len(sys.argv) == 6:
    arg3 = int(sys.argv[4])
    arg4 = int(sys.argv[5])

file_name=str(sys.argv[1])
fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

for i in range(0,N): #only 1 planet for now
    p=data[i::N]
    if arg2 == 3:
        E0 = p[0,arg2]
        y = abs((p[arg4:arg3,arg2] - E0)/E0)
    else:
        y = p[arg4:arg3,arg2]
    plt.plot(p[arg4:arg3,arg1], y, 'o'+colors[i], markeredgecolor='none', ms = 2)#, label='m$_{'+str(i+1)+'}$='+str(round(100*mp[i]/(3*10**(-6)))/100.)+' m$_{earth}$', )

#if arg1==1:
#    plt.ylim([0.699,0.701])
plt.xlim([p[arg4,0],p[arg3,0]])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.show()
