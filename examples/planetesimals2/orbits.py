import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi
import re

names=['time (years)','time (mini, years)','Energy and N_CE', 'r_min','(dt*v_rel/r)_max', 'Energy','Energy and r_min','Energy and (dt*v_rel/r)_max', 'ax','ay','az','par20CE']
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

msval = 2
fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
if arg1 == 6:
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,5], 'o', ms=msval, markeredgecolor='none')
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,3], 'or', ms=msval, markeredgecolor='none')
    plt.xscale('log')
    print 'min_r is', min(data[arg4:arg3,3])
elif arg1 == 7:
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,5], 'o', ms=msval, markeredgecolor='none')
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,4],  'or', ms=msval, markeredgecolor='none')
    plt.xscale('log')
    print 'max(dt*v_rel/r) is', max(data[arg4:arg3,4])
elif arg1 == 2:
    fig, axes = plt.subplots(nrows=2, ncols=1)
    axes[1].plot(data[arg4:arg3,0],data[arg4:arg3,5], 'o', ms=msval, markeredgecolor='none')
    axes[1].set_xscale('log')
    axes[0].plot(data[arg4:arg3,0],data[arg4:arg3,9],  'or', ms=msval, markeredgecolor='none')
    #axes[0].set_xscale('log')
    axes[0].set_ylabel('N_CE')
else:
    plt.plot(data[arg4:arg3,0],data[arg4:arg3,arg1], 'o', ms=msval, markeredgecolor='none')

if arg1 == 5:
    plt.plot(data[arg4:arg3,0],3e-10*data[arg4:arg3,0]**(0.5),color='black',label='t^1/2 growth')
    plt.legend(loc='upper left',prop={'size':10})
plt.xscale('log')

plt.ylabel(names[arg1])
plt.xlabel('time (years)')
plt.yscale('log')
file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
