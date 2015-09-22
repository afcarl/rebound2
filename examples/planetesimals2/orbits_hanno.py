import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names=['time (years)','Semi-Major Axis (AU)','Eccentricity','(Ei - E0) / E0','(Ki - K0) / K0','(Ui - U0) / U0','Total Ang. Mom.','planet-star distance']
colors=['b','g','m','r','c','y']

file_name=str('energy.txt')

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=' ')
plt.plot(data[:,0], data[:,1], 'o', marker='o', markersize=2,markeredgecolor='none',label='energy')
plt.yscale('log')

plt.show()
