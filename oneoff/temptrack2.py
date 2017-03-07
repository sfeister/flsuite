#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
temptrack2.py: Track the temperature, Brms, etc. evolution through time in two simulations

Created by Scott Feister on Wed Dec 14 15:22:10 2016

Follow up with 200kj3.py post analysis

USAGE:
mpirun -np 12 python temptrack2.py

CHANGELOG:
2016-12-14 Created temptrack1.py, based loosely on magtrack3.py
2017-02-21 Created temptrack2.py for updated analysis with RMS
"""

# SHOULD BE FALSE unless this is the most simple re-analysis
quick = False # If true, assume analysis has already been completed and only remake plots

import matplotlib as mpl
mpl.use("Agg")
import numpy as np
import os
import matplotlib.pyplot as plt
from mpi4py import MPI
import cPickle as pickle
import flsuite.sftools as sf
import scipy.constants as sc

if not quick: # Do the full analysis
    import yt
    yt.enable_parallelism() # Tap into yt's mpi4py parallelism (e.g. now can call via mpirun -np 10 python <blah>.py)
    yt.funcs.mylog.setLevel(30) # This sets the output notification threshold to 30, WARNING. Default is 20, INFO.
    from flsuite.flyt import tstools as tst
    from flsuite.flyt import flyt
    import flsuite.flyt.morefields # Add a few extra, custom derived fields to FLASH datasets loaded through YT


def anlzDcentsph(ds, anlsD, rads_um, fld_ids):
    """ Average fields within a spherical ball, of various fixed radii, centered on domain """
    prefix = 'centsph_'
    anlsD['centsph_fld_ids'] = fld_ids # Sphere fields for which to save
    anlsD['centsph_rads_um'] = rads_um # Sphere radii to analyze

    for rad_um in rads_um:
        reg = ds.sphere(ds.domain_center, (rad_um, 'um')) # Centered on domain
        for fld_id in fld_ids:
            anlsD['centsph_mean_' + str(int(rad_um)) + 'um_' + fld_id] = reg.mean(fld_id).v
            anlsD['centsph_rms_' + str(int(rad_um)) + 'um_' + fld_id] =  np.sqrt(np.sum(reg[fld_id]**2)/len(reg[fld_id])).v
            anlsD['centsph_max_' + str(int(rad_um)) + 'um_' + fld_id] =  reg.max(fld_id).v
            anlsD['centsph_std_' + str(int(rad_um)) + 'um_' + fld_id] =  np.std(reg[fld_id]).v
        
    return anlsD

def anlzDcentcube(ds, anlsD, edges_um, fld_ids):
    """ Average fields within a cube, of various fixed edge lengths, centered on domain """
    prefix = 'centcube_'
    anlsD['centcube_fld_ids'] = fld_ids # Cube fields for which to save
    anlsD['centcube_edges_um'] = edges_um # Edge lengths to analyze

    for edge_um in edges_um:
        ledge = ds.domain_center - ds.arr([edge_um/2., edge_um/2., edge_um/2.], 'um') # Cube left edge
        redge = ds.domain_center + ds.arr([edge_um/2., edge_um/2., edge_um/2.], 'um') # Cube right edge
        reg = ds.region(ds.domain_center, ledge, redge) # Cube of edge length edge_um", centered on domain
        for fld_id in fld_ids:
            anlsD['centcube_mean_' + str(int(edge_um)) + 'um_' + fld_id] = reg.mean(fld_id).v
            anlsD['centcube_rms_' + str(int(edge_um)) + 'um_' + fld_id] =  np.sqrt(np.sum(reg[fld_id]**2)/len(reg[fld_id])).v
            anlsD['centcube_max_' + str(int(edge_um)) + 'um_' + fld_id] =  reg.max(fld_id).v
            anlsD['centcube_std_' + str(int(edge_um)) + 'um_' + fld_id] =  np.std(reg[fld_id]).v
        
    return anlsD

def plotTcentsph(anlsT, outdir='.'):
    """ Generically plot the summarized results from anlzDcentsph """
    print("Making central sphere plots...")
    simname = os.path.basename(os.path.normpath(outdir)) # Use the name of the output directory as a label for this sim., to use in plots etc.

    rads_um = anlsT['centsph_rads_um'][0] # Available sphere radii
    fld_ids = anlsT['centsph_fld_ids'][0] # Available field IDs
    tgv = anlsT['times_ns']

    for fld_id in fld_ids:
        fig = plt.figure(1)
        fig.clear()
        ax = fig.add_subplot(111)
        for rad_um in rads_um:
            yvals = anlsT['centsph_mean_' + str(int(rad_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="r=" + str(int(rad_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("Mean of volume: " + fld_id + " (yt units)")
        ax.set_title('Mean of Central sphere (radii r) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central sphere (Mean) - ' + fld_id + '.png'))

        fig = plt.figure(2)
        fig.clear()
        ax = fig.add_subplot(111)
        for rad_um in rads_um:
            yvals = anlsT['centsph_rms_' + str(int(rad_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="r=" + str(int(rad_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("RMS of volume: " + fld_id + " (yt units)")
        ax.set_title('RMS of Central sphere (radii r) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central sphere (RMS) - ' + fld_id + '.png'))

        fig = plt.figure(3)
        fig.clear()
        ax = fig.add_subplot(111)
        for rad_um in rads_um:
            yvals = anlsT['centsph_max_' + str(int(rad_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="r=" + str(int(rad_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("Max in volume: " + fld_id + " (yt units)")
        ax.set_title('Max in Central sphere (radii r) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central sphere (Max) - ' + fld_id + '.png'))

def plotTcentcube(anlsT, outdir='.'):
    """ Generically plot the summarized results from anlzDcentcube """
    print("Making central cube plots...")
    simname = os.path.basename(os.path.normpath(outdir)) # Use the name of the output directory as a label for this sim., to use in plots etc.

    edges_um = anlsT['centcube_edges_um'][0] # Available cube edge lengths
    fld_ids = anlsT['centcube_fld_ids'][0] # Available field IDs
    tgv = anlsT['times_ns']

    for fld_id in fld_ids:
        fig = plt.figure(1)
        fig.clear()
        ax = fig.add_subplot(111)
        for edge_um in edges_um:
            yvals = anlsT['centcube_mean_' + str(int(edge_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="l=" + str(int(edge_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("Mean of volume: " + fld_id + " (yt units)")
        ax.set_title('Mean of Central cube (Edge length l) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central cube (Mean) - ' + fld_id + '.png'))

        fig = plt.figure(2)
        fig.clear()
        ax = fig.add_subplot(111)
        for edge_um in edges_um:
            yvals = anlsT['centcube_rms_' + str(int(edge_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="l=" + str(int(edge_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("RMS of volume: " + fld_id + " (yt units)")
        ax.set_title('RMS of Central cube (Edge length l) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central cube (RMS) - ' + fld_id + '.png'))

        fig = plt.figure(3)
        fig.clear()
        ax = fig.add_subplot(111)
        for edge_um in edges_um:
            yvals = anlsT['centcube_max_' + str(int(edge_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="l=" + str(int(edge_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("Max in volume: " + fld_id + " (yt units)")
        ax.set_title('Max in Central cube (Edge length l) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central cube (Max) - ' + fld_id + '.png'))

        fig = plt.figure(4)
        fig.clear()
        ax = fig.add_subplot(111)
        for edge_um in edges_um:
            yvals = anlsT['centcube_std_' + str(int(edge_um)) + 'um_' + fld_id]
            ax.plot(tgv, yvals, label="l=" + str(int(edge_um)) + " um")
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("StDev in volume: " + fld_id + " (yt units)")
        ax.set_title('StDev in Central cube (Edge length l) vs. Time: ' + fld_id)
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, 'Central cube (StDev) - ' + fld_id + '.png'))
        
def anlzD(ds, outdir='.'):
    """Custom single-step analysis callback"""
    ## Initialize anlsD and add entry for current time
    anlsD = {}
    anlsD['times_ns'] = ds.current_time.in_units('ns').v
    #anlsD = anlzDcentsph(ds, anlsD, [500, 1000, 2000], ['tele', 'magtot', 'magsqr', 'KEdens', 'Pm', 'dens', 'veltot'])
    anlsD = anlzDcentcube(ds, anlsD, [1000, 2000, 4000], ['tele', 'magtot', 'magsqr', 'KEdens', 'Pm', 'dens', 'veltot'])
    
    return anlsD

def plotT(anlsT, outdir='.'):
    """Custom plotting callback"""
    
    #plotTcentsph(anlsT, outdir=outdir)
    plotTcentcube(anlsT, outdir=outdir)
    
    return anlsT

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    #simname = 'NIF_TDYNO_200KJ'
    #datdir = r'/projects/CosmicLaser/tzeferac/SCRIPT1/RUN3'
    #basenm = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
    
    simname = "NIF_TDYNO_BAND"
    datdir = r'/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND'
    basenm = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file

    fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_????') # yt-time-series filename pattern for plot files
    # hdf5_plt_cnt_????
    # hdf5_plt_cnt_???[0,2,4,6,8]
    # hdf5_plt_cnt_???[0,5]
    # hdf5_plt_cnt_??[0,2,4,6,8]0
    # hdf5_plt_cnt_??[0,5]0
    outroot=r'/home/sfeister/myouts/feb2017_NIF_comparison'
    
    if rank == 0:
        print("STARTING WORK ON: " + simname)
        outdir = sf.subdir(outroot, simname)
    else:
        outdir = None
    outdir = comm.bcast(outdir, root=0)
    
    if quick: # Quick and dirty re-analysis
        with open(os.path.join(outdir, "anlsT.p")) as f:
            anlsT = pickle.load(f)
        plotT(anlsT, outdir=outdir)
    else: # Full analysis
        # Subfolder for this particular analysis run
        ts = yt.load(fnpatt)
        anlsT = tst.anlzT(ts, anlzD, outdir=outdir, plotT=plotT)
