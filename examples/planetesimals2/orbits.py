import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names=['time (years)','Semi-Major Axis (AU)','Eccentricity','(Ei - E0) / E0','(Ki - K0) / K0','(Ui - U0) / U0','Total Ang. Mom.','planet-star distance']
colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

#Get number of massive planets
N_active = int(sys.argv[3])

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=0
arg2=int(sys.argv[2])

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

for i in xrange(0,N_active):
    p=data[i::N_active]
    plt.plot(p[:,arg1], p[:,arg2], 'o'+colors[i], marker='o', markersize=2,markeredgecolor='none',label='planet '+str(i))
    if arg2 == 3 or arg2==4 or arg2==5 or arg2==6:
        plt.yscale('log')
        break

#plt.xlim([33,45])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.show()
