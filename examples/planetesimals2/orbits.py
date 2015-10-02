import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

#names=['time (years)','Semi-Major Axis (AU)','Eccentricity','(Ei - E0) / E0','(Ki - K0) / K0','(Ui - U0) / U0','Total Ang. Mom.','planet-star distance']
names=['time (years)','time (mini, years)','N_encounter_previous', 'r_min','(dt*v_rel/r)_max', 'Energy','Kinetic','Potential','Energy and r_min','Energy and (dt*v_rel/r)_max']
colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=int(sys.argv[2])
#arg2=int(sys.argv[3])
arg3 = -1
arg4 = 0
if len(sys.argv) >= 4:
    arg3 = int(sys.argv[3])
if len(sys.argv) >= 5:
    arg4 = int(sys.argv[4])

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
if arg1 == 8:
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,5])
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,3], 'r')
elif arg1 == 9:
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,5])
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,4], 'r')
else:
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,arg1])

plt.ylabel(names[arg1])
plt.xlabel('time (years)')
plt.yscale('log')
plt.show()

#Get number of massive planets
#fos = open(file_name[0:-4]+'_Properties.txt', 'r')
#for i in xrange(0,3):
#    header = fos.readline()
#output = header.split(",")
#N_active = int(output[2]) - 1    #-1 cause of the star
#N_active = 2

#for i in xrange(0,N_active):
#    p=data[i::N_active]
#    plt.plot(p[arg4:arg3,arg1], p[arg4:arg3,arg2], 'o'+colors[i], marker='o', markersize=2,markeredgecolor='none',label='planet '+str(i), )
#    if arg2 == 3 or arg2==4 or arg2==5 or arg2==6:
#        plt.yscale('log')
#        break

#plt.ylim([1e-4,0.05])
#if arg2==7:
#plt.ylim([0.69,0.71])
    #plt.ylim([0.98,1.02])
    #if arg2==1:
#plt.ylim([0.699,0.701])
    #plt.ylim([0.998,1.002])
#plt.xlim([p[arg4,0],p[arg3,0]])
#plt.xlim([33,45])
#plt.xlabel('' + names[arg1])
#plt.ylabel('' + names[arg2])
#plt.legend(loc='upper left',prop={'size':10})
