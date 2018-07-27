#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tsexample2.py: Example, Make 2D yt slice plots at every timestep (a png movie)

Modify the file pattern for your HDF5 files ("filepatt" variable below), then run via "python tsexample2.py".
If mpi4py is installed on your system, you could run this in parallel on four cores by "mpirun -np 4 python tsexample2.py"

Created by Scott Feister on Fri Jul 27 09:38:35 2018
"""

import os, sys
import yt
import flsuite.flyt.tstools as tst
import matplotlib.pyplot as plt

def anlzD(ds, anlsD, outdir):
    """ Use yt to make plots (z-axis slice) of density and electron temperature for this timestep """
    pltnum = str(ds)[-4:] # Plot number string, e.g. strips "0020" off of "lasslab_hdf5_plt_cnt_0020"

    slc = yt.SlicePlot(ds, 'z', 'dens')
    slc.save('dens_' + pltnum + ".png") # E.g. "dens_0020.png"

    slc = yt.SlicePlot(ds, 'z', 'tele')
    slc.save('pres_' + pltnum + ".png") # E.g. "pres_0020.png"
    
if __name__ == "__main__":
    filepatt = '/home/scott/myouts/lasslab_hdf5_plt_cnt_????' # File pattern (note the ????)
    ts = yt.load(filepatt)
    anlsT = tst.anlzT(ts, anlzD, outdir='.')
