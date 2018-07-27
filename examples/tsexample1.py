#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tsexample1.py: Example, Make 1D line plot of maximum density vs. time

Modify the file pattern for your HDF5 files ("filepatt" variable below), then run via "python tsexample1.py".
If mpi4py is installed on your system, you could run this in parallel on four cores by "mpirun -np 4 python tsexample1.py"

Created by Scott Feister on Fri Jul 27 09:38:35 2018
"""

import os, sys
import yt
import flsuite.flyt.tstools as tst
import matplotlib.pyplot as plt

def anlzD(ds, anlsD, outdir):
    """ Save timestamp and max density for this dataset """
    anlsD['time_ns'] = ds.current_time.in_units('ns') # A single value
    anlsD['maxdens'] = ds.find_max('dens')[0] # A single value
    
def plotT(anlsT, outdir):
    """ Plot max density vs. time over all datasets """
    fig, ax = plt.subplots()
    ax.plot(anlsT['time_ns'], anlsT['maxdens']) # A 1D array, and a 1D array
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Density (g/cc)")
    ax.set_title("Maximum density in simulation domain vs. Time")
    fig.savefig(os.path.join(outdir, 'densvtime.png'))
    plt.close(fig)
    
if __name__ == "__main__":
    filepatt = '/home/scott/myouts/lasslab_hdf5_plt_cnt_????' # File pattern (note the ????)
    ts = yt.load(filepatt)
    anlsT = tst.anlzT(ts, anlzD, outdir='.', plotT=plotT)
