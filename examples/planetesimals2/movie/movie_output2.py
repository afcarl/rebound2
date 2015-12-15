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

N_massive = int(raw_input("Number of massive bodies (including sun): ")) - 1

dir = 'movie_output/'
files = glob.glob(dir+'hybridbody*.txt')
files = sorted(files, key=natural_key)
N_bodies = len(files)  #number of files we're dealing with
plotrange = 1
N_prev = 0

#read in data for each body
data = []
try:
    for f in files:
        ff = open(f, 'r')
        lines = ff.readlines()
        data.append(lines)
except:
    print 'couldnt open', f, 'exiting'
    exit(0)

n_it = len(data[0])   #calc num of lines
limit = 1               #size limits for plots = (x,y,z/2)

print 'deleting any existing .png images in output_movie folder'
call("rm output_movie/*.png",shell=True)

collision = np.zeros(N_bodies)
for iteration in xrange(0,n_it):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in xrange(0,N_bodies):
        try:
            line = data[i][iteration].split(",")
            color = get_color(int(line[1]),N_massive)
            ax.scatter(float(line[2]),float(line[3]),float(line[4]),c=color, lw=0)
        except:
            if collision[i] == 0:
                print 'particle',i,' had a collision'
                collision[i] = 1
    #plotting details - make it all look pwetty.
    ax.set_xlim([-limit,limit])
    ax.set_ylim([-limit,limit])
    ax.set_zlim([-limit/4,limit/4])
    ax.view_init(elev = 90, azim=100)    #both are in degrees. elev = 0 or 90 is what you want
    output = 't='+line[0]
    ax.set_title(output)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.savefig('movie_output/movie_output'+str(iteration)+'.png')
    print 'completed iteration '+str(iteration + 1)+' of '+str(n_it)
