#This macro plots the distribution of close encounters for a given run. Takes *_CE.txt as input

import sys
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import re

filename=str(sys.argv[1])
f = open(filename, 'r')
lines = f.readlines()
CE_times = np.zeros(0)

title = filename.split("_")
tit_length = len(title)

for line in lines:
    if line.startswith('t=',0,2):
        split = line.split(',')
        CE_times = np.append(CE_times, float(split[0][2:]))

pl.hist(CE_times, bins = np.logspace(-1,np.log10(float(title[tit_length-4][1:])), num=100), color='red')
pl.gca().set_xscale("log")
pl.gca().set_xlabel("time (years)")
pl.gca().set_ylabel("Number of Close Encounters")
pl.gca().set_title("Distribution of Close Encounters for "+title[tit_length-3]+", "+title[tit_length-2])
file_output_name = re.sub('\.txt$', '', filename)
plt.savefig(file_output_name+'.png')
pl.show()