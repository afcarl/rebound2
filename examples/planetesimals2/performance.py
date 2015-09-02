#The purpose of this macro is to compare the performance of my Hybrid integrator to the IAS15

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
pi = math.pi

dir = 'saved_output/t=5000/'

IASfile = 'planet_IAS15_t=5000.txt'
HYBfile = 'planet_HYBRID_5_1.txt'

names=['time (years)','Semi-Major Axis (AU)','Eccentricity','(Ei - E0) / E0','Total Ang. Mom.','planet-star distance']
colorsHYB=['b','g','m','r','c','y']
colorsIAS='black'

fos = open(dir+HYBfile[0:-4]+'_Properties.txt', 'r')
for i in xrange(0,3):
    header = fos.readline()
output = header.split(",")
N_active = int(output[2]) - 1    #-1 cause of the star

arg1=int(sys.argv[1])
arg2=int(sys.argv[2])
arg3 = 'both'
if len(sys.argv) == 4:
    arg3 = sys.argv[3]

#open files
fosHYB = open(dir+HYBfile, 'r')
dataHYB = np.loadtxt(fosHYB, delimiter=',')
fosIAS = open(dir+IASfile, 'r')
dataIAS = np.loadtxt(fosIAS, delimiter=',')

lo = 0
hi = N_active
if arg3 == 'inner':
    lo = 0
    hi = 1
if arg3 == 'outer':
    lo = N_active - 1
    hi = N_active

for i in xrange(lo,hi):
    p=dataHYB[i::N_active]
    q=dataIAS[i::N_active]
    plt.plot(q[:,arg1], q[:,arg2], 'ok', markeredgecolor='none', ms = 3, label='IAS planet '+str(i+1), )
    plt.plot(p[:,arg1], p[:,arg2], 'o'+colorsHYB[i], markeredgecolor='none', ms = 2, label='HYBRID planet '+str(i+1), )
    if arg2 == 3 or arg2==4:
        plt.yscale('log')
        break

if arg2==5:
    if arg3 == 'inner':
        plt.ylim([0.69,0.71])
    if arg3 == 'outer':
        plt.ylim([0.99,1.01])
if arg2==1:
    if arg3 == 'inner':
        plt.ylim([0.699,0.701])
    if arg3 == 'outer':
        plt.ylim([0.999,1.001])
plt.xlim([0,p[-1,0]])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.show()