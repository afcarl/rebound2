#This macro plots the distribution of ejections/collisions for a given run(s). Can also compare identical runs from Swifter and my Hybrid integrator. Takes *_removed.txt as input, or the corresponding file for swifter.

import sys
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import re

swifter_compare = 0

filename=str(sys.argv[1])
f = open(filename, 'r')
lines = f.readlines()
eject_times = np.zeros(0)
collide_times = np.zeros(0)

title = filename.split("_")
tit_length = len(title)
skipped_lines = 0

for line in lines:
    try:
        split = line.split()
        if split[0] == "Collision":
            try:
                collide_times = np.append(collide_times, float(re.sub('t=','',split[2])))
            except:
                collide_times = np.append(collide_times, float(split[3]))
        elif split[0] == "Ejection":
            try:
                eject_times = np.append(eject_times, float(re.sub('t=','',split[2])))
            except:
                eject_times = np.append(eject_times, float(split[3]))
    except:
            skipped_lines += 1

try:
    outerlim = float(title[tit_length-4][1:])
    title1 = title[tit_length-3]
    title2 = title[tit_length-2]
except:
    outerlim = float(raw_input("Enter length of simulation: "))
    title1 = raw_input("Enter Title (e.g. Symba): ")
    title2 = raw_input("Enter Number of Planetesimals (e.g. Np500): ")

pl.hist(eject_times, bins = np.logspace(-1,np.log10(outerlim), num=50), color='orange', alpha=0.7,label='Ejected particles, N='+str(len(eject_times)))
pl.hist(collide_times, bins = np.logspace(-1,np.log10(outerlim), num=50), color='green', alpha=0.7, label='Collided particles, N='+str(len(collide_times)))

if swifter_compare == 1:
    try:
        filename2 = str(sys.argv[2])
    except:
        print 'You must include location of swifter file'
        exit(0)
    f = open(filename2, 'r')
    lines = f.readlines()
    eject_times = np.zeros(0)
    collide_times = np.zeros(0)
    for line in lines:
        try:
            split = line.split()
            if split[0] == "Collision":
                collide_times = np.append(collide_times, float(split[3]))
            elif split[0] == "Ejection":
                eject_times = np.append(eject_times, float(split[3]))
        except:
            skipped_lines += 1
    pl.hist(eject_times, bins = np.logspace(-1,np.log10(outerlim), num=50), color='red', alpha=0.5, label='Swifter Ejected particles, N='+str(len(eject_times)))
    pl.hist(collide_times, bins = np.logspace(-1,np.log10(outerlim), num=50), color='blue', alpha = 0.75, label='Swifter Collided particles, N='+str(len(collide_times)))

pl.legend(loc='upper left',prop={'size':10})
pl.gca().set_xscale("log")
pl.gca().set_xlabel("time (years)")
pl.gca().set_ylabel("Number of Ejected Planetesimals")
pl.gca().set_title("Distribution of Planetesimal Ejections for "+title1+", "+title2)
file_output_name = re.sub('\.txt$', '', filename)
plt.savefig(file_output_name+'_Removedplot.png')
pl.show()