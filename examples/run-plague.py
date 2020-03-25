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
    os.system("Rscript plague-example.R year 2017 nMCMC 2e6 epi_model 2 Temp 10 "+" nfit "+str(nfit))

nthreads = 48
nfit = np.arange(57,105,dtype=np.integer)


if __name__ == '__main__':
    info('main line')
    p = multiprocessing.Pool(nthreads)
    p.map(f, nfit)

