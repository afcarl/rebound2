#This macro plots the distribution of ejections for a given run. Takes *.txt as input

import sys
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import re

filename=str(sys.argv[1])
f = open(filename, 'r')
lines = f.readlines()
eject_times = np.zeros(0)

title = filename.split("_")
tit_length = len(title)

split = lines[0].split(",")
N = int(split[2])
for line in lines:
    split = line.split(",")
    if int(split[2]) < N:
        eject_times = np.append(eject_times, float(split[0]))
        N = int(split[2])

pl.hist(eject_times, bins = np.logspace(-1,np.log10(float(title[tit_length-3][1:])), num=50), color='orange')
pl.gca().set_xscale("log")
pl.gca().set_xlabel("time (years)")
pl.gca().set_ylabel("Number of Ejected Planetesimals")
pl.gca().set_title("Distribution of Planetesimal Ejections for "+title[tit_length-2]+", "+title[tit_length-1])
file_output_name = re.sub('\.txt$', '', filename)
plt.savefig(file_output_name+'_eject.png')
pl.show()