#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
#import numpy as np

#params=[(50000,10,10),(50000,50,11),(50000,100,12),(50000,500,13),(50000,1000,14),(50000,5000,15),(1000,100,16),(5000,100,17),(10000,100,18),(100000,100,19),(500000,100,21)]
#params=[(5000000,100,12),(5000000,1000,21),(100000,10,30),(100000,100,13),(100000,1000,14),(100000,5000,54),(100000,10000,80)]
params=[(3000000,100,1e-5,94),(3000000,100,1e-7,94),(3000000,100,1e-9,94),(1000000,1000,1e-5,12),(1000000,1000,1e-6,12),(1000000,1000,1e-7,12),(1000000,1000,1e-9,12),(100000,10000,1e-7,44),(100000,10000,1e-9,44)]
length = len(params)

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[params[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()