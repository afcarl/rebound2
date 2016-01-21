#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
#import numpy as np

#params=[(10000000,50,13),(10000000,50,26),(10000000,50,34),(10000000,50,51),(10000000,50,63),(10000000,50,65),(10000000,50,70),(10000000,50,10),(5000000,100,16),(5000000,100,18),(5000000,100,28),(5000000,100,37),(5000000,100,52),(5000000,100,62),(5000000,100,67),(5000000,100,72),(1000000,500,11),(1000000,500,12),(1000000,500,22),(1000000,500,33),(1000000,500,49),(1000000,500,62),(1000000,500,65),(1000000,500,71)]
params=[(5000000,100,16),(5000000,100,18),(5000000,100,28),(5000000,100,37),(5000000,100,52),(5000000,100,62),(5000000,100,67),(5000000,100,72)]
length = len(params)

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[params[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
