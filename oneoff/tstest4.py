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
    slc = yt.SlicePlot(ds, 'z', ['dens'])
    slc.save(outdir)
    
    anlsD = {}
    anlsD['times_ns'] = ds.current_time.in_units('ns').v

    xgv, densgv = flyt.lineout(ds, field='dens', axis=0)
#    print ds.domain_center[1]
#    print ds.domain_center[2]
#    print type(ds.domain_center[2])
#    ray = ds.ortho_ray(0, (0.5, 0.5))
#    srt = np.argsort(ray['x'])
#    rayx = np.array(ray["x"][srt])
#    rayfld = np.array(ray["dens"][srt])
#    
#    #dx = np.max(np.diff(rayx))/4.0 # Somewhat arbitrary
#    xgv = np.linspace(np.min(rayx), np.max(rayx), 200)
#    fldgv = np.interp(xgv, rayx, rayfld)
    anlsD['x_line'] = xgv
    anlsD['dens_line'] = densgv
    
    return anlsD

def plotT(anlsT, outdir='.'):
    plt.figure(1)
    plt.clf()
    plt.plot(anlsT['x_line'].T, anlsT['dens_line'].T)
    plt.savefig(os.path.join(outdir, "densplot.png"))
    
    plt.figure(2)
    plt.clf()
    xgv = anlsT['x_line'][0]*1e4 # xgv in microns
    tgv = anlsT['times_ns']
    ax = plt.subplot(111)
    ax.pcolorfast((tgv[0], tgv[-1]), (xgv[0], xgv[-1]), anlsT['dens_line'].T)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("X (um)")
    plt.savefig(os.path.join(outdir, "timeplot.png"))

    return anlsT, xgv, tgv

if __name__ == "__main__":
    fnpatt = r"C:\Users\Scott\Documents\temp\TDYNO vids\windtunnel test\windtunnel_4lev_hdf5_plt_cnt_????"
    outdir = r"C:\Users\Scott\Documents\temp\TDYNO vids\windtunnel test\outs"
    #ts = yt.load(fnpatt)
    #anlsT = tst.anlzT(ts, anlzD, outdir=outdir, plotT=plotT)
