import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import matplotlib.cm as cm
pi = math.pi
from mpl_toolkits.mplot3d import Axes3D
import re

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

file_names = 'xyz_temp/outOct28mini_'
files = glob.glob(file_names+'*.txt')
files = sorted(files, key=natural_key)
N_files = len(files)  #number of files we're dealing with
N_planets = 500

#colors
colors =  ["black" for x in range(N_planets)]
colors[0] = 'yellow'
colors[1] = 'orange'
colors[2] = 'red'

for i in xrange(0,N_files):
    fos = open(files[i], 'r')
    t,n,x,y,z = np.loadtxt(fos, delimiter=',', unpack='True')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,c=colors[0:len(x)], lw=0)
    ax.scatter(x[0:3],y[0:3],z[0:3],c=colors[0:3],lw=0, s=30)
    #ax.set_title('t='+str(t[0])+', t$_{error}$ = 201.826486')
    ax.set_title('t='+str(t[0])+', t$_{error}$ = 48.297440')
    plt.savefig(file_names+''+str(i)+'.pdf')
    print 'completed iteration',i
