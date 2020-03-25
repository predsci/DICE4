#!/usr/bin/python

import multiprocessing
from multiprocessing import Process
import subprocess
import os
import string
import random
import numpy as np

def info(title):
    print title
    print 'module name:',__name__
    if hasattr(os,'getppid'): 
         print 'parent process:',os.getpid()
    print 'process id: ',os.getpid()

def f(nfit):
    os.system("Rscript sd-example.R year 2017 nMCMC 5e6 epi_model 2 model 5 isingle 1 Temp 2 "+" nfit "+str(nfit))

nthreads = 4
nfit = np.arange(40,44,dtype=np.integer)

#nthreads = 4
#nfit = np.arange(49,53,dtype=np.integer)

if __name__ == '__main__':
    info('main line')
    p = multiprocessing.Pool(nthreads)
    p.map(f, nfit)

