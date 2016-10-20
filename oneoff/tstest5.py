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
import flsuite.sftools as sf
from flsuite.flyt import flyt
#import flsuite.flyt.morefields # Add a few extra, custom derived fields to FLASH datasets loaded through YT
import numpy as np
import os
import matplotlib.pyplot as plt

def anlzD(ds, outdir='.'):
#    slc = yt.SlicePlot(ds, 'z', ['velz'])
#    slc.save(outdir)
#    slc = yt.SlicePlot(ds, 'x', ['velz'])
#    slc.save(outdir)
#    slc = yt.SlicePlot(ds, 'y', ['velz'])
#    slc.save(outdir)

    anlsD = {}
    anlsD['times_ns'] = ds.current_time.in_units('ns').v
    anlsD['z_line'], anlsD['velz_line'] = flyt.lineout(ds, field='velz', axis=2, npts=2000)

    return anlsD

def plotT(anlsT, outdir='.'):
  
    fig = plt.figure(2)
    plt.clf()
    xgv = anlsT['z_line'][0]*1e1 # xgv in mm
    tgv = anlsT['times_ns']
    vmax = 2*np.mean(np.abs(anlsT['velz_line']))
    vmin = -vmax
    ax = plt.subplot(111)
    cax = ax.pcolorfast((tgv[0], tgv[-1]), (xgv[0], xgv[-1]), anlsT['velz_line'].T, vmin=vmin, vmax=vmax, cmap='RdBu')
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Z (mm)")
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel("Velz (linear, a.u.)")
    ax.set_title("Velz lineout X=Y=0, through time")
    plt.savefig(os.path.join(outdir, "timeplot.png"), dpi=300)

    return anlsT

if __name__ == "__main__":
    # Note that this will run in parallel, so no prints please
    datdir = r'/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND' # Directory holding HDF5 FLASH outputs
    basenm = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
    fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_????') # yt-time-series filename pattern for plot files

    outdir = r'/home/sfeister/myouts/TDYNO_BAND' # Folder for outputs
    ts = yt.load(fnpatt)
    anlsT = tst.anlzT(ts, anlzD, outdir=outdir, plotT=plotT)

