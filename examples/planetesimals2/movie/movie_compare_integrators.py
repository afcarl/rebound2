#This integrator will compare the hybrid integrator output to the swifter integrator output. Basic instructions: 1) Run a simulation with your hybrid integrator (no movie). 2) Run the comparable one with swifter. 3) Check the output of swifter to get the dt increment, so that you can set that dt increment to be your movie output frequency. Then put the output files in the movie_compare_integrators folder and run this.

import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import math
pi = math.pi
from mpl_toolkits.mplot3d import Axes3D
import re
from subprocess import call

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def get_color(id,N_massive):
    color = 'black'
    if id == 0:
        color = 'yellow'
    elif id <= N_massive:
        color = 'red'
    return color

dir = 'movie_compare_integrators/'
fileshyb = glob.glob(dir+'hybridbody*.txt')
fileshyb = sorted(fileshyb, key=natural_key)
files2 = glob.glob(dir+'follow*.txt')
files2 = sorted(files2, key=natural_key)
N_bodies = len(files2)  #number of files we're dealing with

#read in data for each body - hyb
datahyb = []
try:
    for f in fileshyb:
        ff = open(f, 'r')
        lines = ff.readlines()
        datahyb.append(lines)
except:
    print 'couldnt open', f, 'exiting'
    exit(0)

#read in data for each body - the second integrator in question
data2 = []
try:
    for f in files2:
        ff = open(f, 'r')
        lines = ff.readlines()
        data2.append(lines)
except:
    print 'couldnt open', f, 'exiting'
    exit(0)

n_it = len(datahyb[0])     #calc num of lines

print 'deleting any existing .png images in output_movie folder'
call("rm output_movie/*.png",shell=True)

collision = np.zeros(N_bodies)

xyz_space = int(sys.argv[1])

#dx,dy,dz vs. time for each set of objects
if xyz_space == 0:
    for i in xrange(1,N_bodies):
        dx = np.zeros(n_it)
        dy = np.zeros(n_it)
        dz = np.zeros(n_it)
        t = np.zeros(n_it)
        for j in xrange(0,n_it):
            if collision[i] == 0:
                try:
                    linehyb = datahyb[i][j].split(",")
                    line2 = data2[i][j].split()
                    t[j] = float(linehyb[0])
                    dx[j] = abs(float(line2[2])-float(linehyb[2]))
                    dy[j] = abs(float(line2[3])-float(linehyb[3]))
                    dz[j] = abs(float(line2[4])-float(linehyb[4]))
                except:
                    if collision[i] == 0:
                        print 'particle',i,' had a collision'
                        collision[i] = 1
        plt.plot(t, dx, 'o', color='blue',label='dx')
        plt.plot(t, dy, 'o', color='green',label='dy')
        plt.plot(t, dz, 'o', color='red',label='dz')
        plt.legend(loc='upper left',prop={'size':10})
        plt.yscale('log')
        plt.xlabel('time (years)')
        plt.ylabel('difference between hybrid and competitor integration')
        plt.savefig('movie_compare_integrators/MCO'+str(i)+'.png')
        plt.clf()
        print 'finished body', i

#This is the plotting in x,y,z space
else:
    limit = 1e-6               #size limits for plots = (x,y,z/2)
    coloring = ['yellow','red','green','blue','purple']
    names = ['Sun','Jupiter','Saturn','Uranus','Neptune']
    N_massive = int(raw_input("Number of massive bodies (including sun): ")) - 1
    for iteration in xrange(0,n_it):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(0,0,0,c='yellow', lw=0)  #sun
        for i in xrange(1,N_bodies):
            if collision[i] == 0:
                try:
                    linehyb = datahyb[i][iteration].split(",")
                    line2 = data2[i][iteration].split()
                    color = get_color(int(linehyb[1]),N_massive)
                    ax.scatter(float(line2[2])-float(linehyb[2]),float(line2[3])-float(linehyb[3]),float(line2[4])-float(linehyb[4]),c=coloring[i], lw=0, label=names[i])
                except:
                    if collision[i] == 0:
                        print 'particle',i,' had a collision'
                        collision[i] = 1
        #plotting details - make it all look pwetty.
        ax.set_xlim([-limit,limit])
        ax.set_ylim([-limit,limit])
        ax.set_zlim([-limit/4,limit/4])
        ax.view_init(elev = 90, azim=100)    #both are in degrees. elev = 0 or 90 is what you want
        output = 't='+linehyb[0]
        ax.set_title(output)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.legend(loc='upper left',prop={'size':10})
        plt.savefig('movie_compare_integrators/MCO'+str(iteration)+'.png')
        print 'completed iteration '+str(iteration + 1)+' of '+str(n_it)
