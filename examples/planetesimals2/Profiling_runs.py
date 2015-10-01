#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
#import numpy as np

#params=[(50000,10,10),(50000,50,11),(50000,100,12),(50000,500,13),(50000,1000,14),(50000,5000,15),(1000,100,16),(5000,100,17),(10000,100,18),(100000,100,19),(500000,100,21)]
#params=[(1000000,10,10),(500000,50,20),(100000,100,30),(75000,500,40),(50000,1000,50),(40000,5000,60)]
#params=[(1000,200,12),(1000,201,22),(1000,199,32),(1000,198,42),(1000,202,52),(1000,203,62)]
params=[(1000,100,13),(1000,101,24),(1000,99,35),(1000,98,46),(1000,102,57),(1000,103,68)]
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