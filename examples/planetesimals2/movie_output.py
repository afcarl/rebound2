import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import math
pi = math.pi
from mpl_toolkits.mplot3d import Axes3D
import re

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def get_colors(N_bodies,N_massive):
    colors =  ["black" for x in range(N_bodies+1)]
    colors[0] = 'yellow'
    for i in xrange(1,N_massive):
        colors[i] = 'red'
    return colors

N_massive = int(raw_input("Number of massive bodies (including sun): "))

dir = 'movie_output/'
files = glob.glob(dir+'movie_output*.txt')
files = sorted(files, key=natural_key)
N_files = len(files)  #number of files we're dealing with
plotrange = 1
N_prev = 0

for i in xrange(0,N_files):
    fos = open(files[i], 'r')
    t,n,x,y,z = np.loadtxt(fos, delimiter=',', unpack='True')
    N = len(x)
    if N != N_prev:
        colors = get_colors(N,N_massive)
        N_prev = N
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,c=colors[0:len(x)], lw=0)
    ax.view_init(elev = 90, azim=100)    #both are in degrees. elev = 0 or 90 is what you want
    ax.set_xlim([-plotrange,plotrange])
    ax.set_ylim([-plotrange,plotrange])
    ax.set_zlim([-plotrange/4,plotrange/4])
    ax.set_title('t='+str(t[0]))
    output = files[i]
    output = re.sub('\.txt$', '', output)
    plt.savefig(output+'.png')
    print 'completed iteration',i
