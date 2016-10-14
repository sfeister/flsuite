#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
tstest4.py: Test tstools parallelization

Created by Scott Feister on Fri Oct 14 09:40:27 2016
"""

import yt
yt.enable_parallelism() # Tap into yt's mpi4py parallelism (e.g. now can call via mpirun -np 10 python <blah>.py)
yt.funcs.mylog.setLevel(30) # This sets the output notification threshold to 30, WARNING. Default is 20, INFO.
from flsuite.flyt import tstools as tst

if __name__ == "__main__":
    fnpatt = r"C:\Users\Scott\Documents\temp\TDYNO vids\windtunnel test\windtunnel_4lev_hdf5_plt_cnt_????"
    ts = yt.load(fnpatt)
    pass
