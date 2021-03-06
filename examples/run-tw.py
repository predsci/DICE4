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
    os.system("Rscript taiwan-example.R year 2002 nMCMC 2e6 epi_model 2 isingle 1 "+" nfit "+str(nfit))

nthreads = 28
nfit = np.arange(25,53,dtype=np.integer)

#nthreads = 4
#nfit = np.arange(49,53,dtype=np.integer)

if __name__ == '__main__':
    info('main line')
    p = multiprocessing.Pool(nthreads)
    p.map(f, nfit)

