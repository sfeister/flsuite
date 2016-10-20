#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
magtrack1.py: Analyze magnetic fields, tracking the zero-velocity front

Created by Scott Feister on Fri Oct 19 2016

USAGE:
mpirun -np 12 python magtrack1.py

CHANGELOG:
2016-10-19 Spun off from "tstest5.py", which was to test parallel processing

TODO:
* Organize the calls for lineout analysis of various widths
* Clean up fatline function (from flyt) to reduce memory overhead of multiple calls
* Incorporate a kinetic energy vs. magnetic energy analysis
* Better incorporate simulation name? Not sure.
* Clean up the lineouts plotting outputs
* Update magnetic field energy so it's not vacuum energy (mu_0) but the real mu?
"""


# SHOULD BE FALSE unless this is the most simple re-analysis
quick = True # If true, assume analysis has already been completed and only remake plots

import matplotlib as mpl
mpl.use("Agg")

import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime
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

def anlzD(ds, outdir='.'):
#    slc = yt.SlicePlot(ds, 'z', ['velz'])
#    slc.save(outdir)
#    slc = yt.SlicePlot(ds, 'x', ['velz'])
#    slc.save(outdir)
#    slc = yt.SlicePlot(ds, 'y', ['velz'])
#    slc.save(outdir)

    ## Initialize anlsD and add entry for current time
    anlsD = {}
    anlsD['times_ns'] = ds.current_time.in_units('ns').v
    
    # Extract a lineout at X=Y=0 for velocity
    anlsD['line1_zgv'], anlsD['line1_velz'] = flyt.lineout(ds, field='velz', axis=2, npts=2000)
    _, anlsD['line1_magtot'] = flyt.lineout(ds, field='magtot', axis=2, npts=2000)
    _, anlsD['line1_dens'] = flyt.lineout(ds, field='dens', axis=2, npts=2000)
    _, anlsD['line1_KEdens'] = flyt.lineout(ds, field='KEdens', axis=2, npts=2000)
    _, anlsD['line1_magsqr'] = flyt.lineout(ds, field='magsqr', axis=2, npts=2000)

    dd = ds.all_data()
    anlsD['line1_velz_units'] = str(dd['velz'].units) # TODO: Better. Include the units in the lineout call
    anlsD['line1_zgv_units'] = str(dd['z'].units) #TODO: Better
    anlsD['line1_dens_units'] = str(dd['dens'].units) #TODO: Better
    
    
    ## TODO: Organize this fatter and fatter lineout section. Also, fix the fatline function so it takes all the fields at once (minimize overhead)
    # Extract a fatter lineout at X=Y=0 for similar quantities
    rad_mm = 0.5 # Radius of the lineout cylinder
    anlsD['line_500um_zgv'], anlsD['line_500um_velz'] = flyt.fatline(ds, rad_mm, fld='velz', npts=2000)
    _, anlsD['line_500um_magtot'] = flyt.fatline(ds, rad_mm, fld='magtot', npts=2000)
    _, anlsD['line_500um_dens'] = flyt.fatline(ds, rad_mm, fld='dens', npts=2000)
    _, anlsD['line_500um_KEdens'] = flyt.fatline(ds, rad_mm, fld='KEdens', npts=2000)
    _, anlsD['line_500um_magsqr'] = flyt.fatline(ds, rad_mm, fld='magsqr', npts=2000)

    # Extract an even fatter lineout at X=Y=0 for similar quantities
    rad_mm = 2.0 # Radius of the lineout cylinder
    anlsD['line_2000um_zgv'], anlsD['line_2000um_velz'] = flyt.fatline(ds, rad_mm, fld='velz', npts=2000)
    _, anlsD['line_2000um_magtot'] = flyt.fatline(ds, rad_mm, fld='magtot', npts=2000)
    _, anlsD['line_2000um_dens'] = flyt.fatline(ds, rad_mm, fld='dens', npts=2000)
    _, anlsD['line_2000um_KEdens'] = flyt.fatline(ds, rad_mm, fld='KEdens', npts=2000)
    _, anlsD['line_2000um_magsqr'] = flyt.fatline(ds, rad_mm, fld='magsqr', npts=2000)

    ## Analyze a lineout of choice for a zero-crossing (searching for first in ROI from +z to -z)
    zgv = anlsD['line1_zgv'] # Extract velocity, in cm
    velz = anlsD['line1_velz']
    ct = (zgv > 2.2*1e-1) & (zgv < 5.8*1e-1) # Define ROI: Z limits in which to look for stagnation
    
    mask_cross = np.append(np.diff(np.signbit(velz)), False) # Boolean array of zero crossing (True if crossing occurs at index)
    cross_tupl = np.where(mask_cross & ct)
    if np.sum(mask_cross & ct) > 0: # If any crossings were found in ROI 
        zstag = zgv[mask_cross & ct][-1]
    else: # No crossings found; just use the maximum Z value in the ROI
        zstag = zgv[ct][-1]
    
    anlsD['z_stag'] = zstag # Define the max-z zero-crossing as stagnation z position
    
    ## Perform the magnetic field sphere-averaging analysis; both on the central sphere and outer spheres
    rads_um = [250, 500, 1000, 2000, 4000] # Sphere radii to explore
    anlsD['sph_rads_um'] = np.array(rads_um)
    for rad_um in rads_um:
        reg = ds.sphere(ds.domain_center, (rad_um, 'um')) # Centered on domain
        anlsD['centsph_' + str(int(rad_um)) + 'um_magtot'] = reg.mean("magtot").v
        anlsD['centsph_' + str(int(rad_um)) + 'um_KEdens'] = reg.mean("KEdens").v
        anlsD['centsph_' + str(int(rad_um)) + 'um_magsqr'] = reg.mean("magsqr").v
        anlsD['centsph_' + str(int(rad_um)) + 'um_voltot'] = reg.sum("cell_volume").v
        
        cent = ds.arr([0., 0., zstag], 'code_length') # Centered at (0, 0, zstag)
        reg = ds.sphere(cent, (rad_um, 'um'))
        anlsD['tracksph_' + str(int(rad_um)) + 'um_magtot'] = reg.mean("magtot").v
        anlsD['tracksph_' + str(int(rad_um)) + 'um_KEdens'] = reg.mean("KEdens").v
        anlsD['tracksph_' + str(int(rad_um)) + 'um_magsqr'] = reg.mean("magsqr").v
        anlsD['tracksph_' + str(int(rad_um)) + 'um_voltot'] = reg.sum("cell_volume").v

    return anlsD

def plotTsph(anlsT, outdir='.'):
    """From anlsT, generate spheres plots (field value averaged over various spheres, vs. time)"""
    print("Making sphere plots...")
    simname = os.path.basename(os.path.normpath(outdir)) # Use the name of the output directory as a label for this sim., to use in plots etc.

    tgv = anlsT['times_ns']
    
    ## Magnetic field vs. time plots (traveling and fixed sphere)
    prefixes = ['centsph_', 'tracksph_']
    titles = ["Averaged |B| over spherical volume\n(Fixed at domain center)",
        "Averaged |B| over spherical volume\n (Z origin changing to follow stagnation for X=Y=0)"]
    outnames = ["B vs time - Centered sphere", "B vs time - Traveling sphere"]

    rad_ums = anlsT['sph_rads_um'][0] # The radius of spheres
    for i in range(len(prefixes)):
        fig = plt.figure(1)
        fig.clear()
        ax = fig.add_subplot(111)
        for rad_um in rad_ums:
            k = prefixes[i] + str(int(rad_um)) + 'um_magtot' # Key in anlsT for this particular combo
            B_kG = anlsT[k]*np.sqrt(4*np.pi)*1e-3 # Convert to kG
            ax.plot(tgv, B_kG, label="Sphere, r=" + str(int(rad_um)) + " um")
            #ax.set_xlim(0, 50)
        ax.set_xlabel("Simulation time (ns)")
        ax.set_ylabel("Volume-averaged magnetic field strength (kG)")
        ax.set_title(titles[i])
        plt.legend(loc='upper left') # Object-oriented cop-out!
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, outnames[i] + '.png'))

    ## Magnetic energy vs. time plot (fixed sphere)
    fig = plt.figure(2)
    fig.clear()
    ax = fig.add_subplot(111)
    for rad_um in rad_ums:
        k = 'centsph_' + str(int(rad_um)) + 'um_magsqr' # Key in anlsT for this particular combo
        Bsq_Tsq = (anlsT[k]*4*np.pi)*1e-8 # Convert to Tesla^2 from code_magnetic^2
        sphvol_m3 = anlsT['centsph_' + str(int(rad_um)) + 'um_voltot'] * 1e-6 # Volume that was averaged over, in m^3
        Etot_J = 0.5 * (Bsq_Tsq / sc.mu_0) * sphvol_m3  # Energy in magnetic fields: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/engfie.html
        ax.plot(tgv, Etot_J, label="Sph., r=" + str(int(rad_um)) + " um")
        #ax.set_xlim(0, 50)
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel("Magnetic Energy (J)")
    ax.set_title("Total magnetic field energy through time,\nwithin spherical volume")
    ax.legend(loc='upper left')
    fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'Magnetic E vs time (lin).png'))

    ax.set_ylim(np.array(ax.get_ylim())/20.0)
    fig.savefig(os.path.join(outdir, 'Magnetic E vs time (lin2).png'))
    
    ax.set_yscale('log')
    ax.set_ylim(1e-7, 1e3)
    ax.legend(loc='lower right')
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'Magnetic E vs time (log).png'))

    ## Kinetic vs. time plot (fixed sphere)
    fig = plt.figure(4)
    fig.clear()
    ax = fig.add_subplot(111)
    for rad_um in rad_ums:
        KEdens_ergpcm3 = anlsT['centsph_' + str(int(rad_um)) + 'um_KEdens'] # Kinetic energy density, in erg/cm^3
        sphvol_cm3 = anlsT['centsph_' + str(int(rad_um)) + 'um_voltot'] # Volume that was averaged over, in cm^3
        KEtot_kJ = KEdens_ergpcm3 * sphvol_cm3 * 1e-7 * 1e-3 # Total kinetic energy, in kJ
        ax.plot(tgv, KEtot_kJ, label="Sph., r=" + str(int(rad_um)) + " um")
        #ax.set_xlim(0, 50)
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel("Kinetic Energy (kJ)")
    ax.set_title("Total kinetic energy through time,\nwithin spherical volume")
    ax.legend(loc='upper left')
    fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'KE vs time (lin).png'))

    ax.set_ylim(np.array(ax.get_ylim())/20.0)
    fig.savefig(os.path.join(outdir, 'KE vs time (lin2).png'))
    
    ax.set_yscale('log')
    ax.set_ylim(1e-7, 1e3)
    ax.legend(loc='lower right')
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'KE vs time (log).png'))
    
    ## Kinetic AND Magnetic energy vs. time plot (fixed sphere)
    fig = plt.figure(3)
    fig.clear()
    ax = fig.add_subplot(111)
    for rad_um in rad_ums:
        prfx = 'centsph_' + str(int(rad_um)) + 'um_' # Beginning of the key
        sphvol_cm3 = anlsT[prfx + 'voltot'] # Volume that was averaged over, in cm^3
        Bsq_Tsq = (anlsT[prfx + 'magsqr']*4*np.pi)*1e-8 # Convert to Tesla^2 from code_magnetic^2
        KEdens_ergpcm3 = anlsT[prfx + 'KEdens'] # Kinetic energy density, in erg/cm^3
        
        KEtot_J = KEdens_ergpcm3 * sphvol_cm3 * 1e-7 # Total kinetic energy, in J
        sphvol_m3 = sphvol_cm3 * 1e-6 # Volume that was averaged over, in m^3
        EBtot_J = 0.5 * (Bsq_Tsq / sc.mu_0) * sphvol_m3  # Energy in magnetic fields: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/engfie.html
        ax.plot(tgv, KEtot_J, label="Kinetic, Sph., r=" + str(int(rad_um)) + " um", color='black', linestyle='--')
        ax.plot(tgv, EBtot_J, label="Magnetic, Sph., r=" + str(int(rad_um)) + " um", color='maroon', linestyle='-')
        #ax.set_xlim(0, 50)
    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel("Energy (J)")
    ax.set_yscale('log')
    ax.set_ylim(1e-7, 1e6)
    ax.xaxis.grid(True)
    ax.set_title("Kinetic (dashed), Magnetic (solid) energy through time,\nwithin spherical volume")
    #ax.legend(loc='upper right', bbox_to_anchor=(1, 0.5))
    fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'Energy vs time - Kinetic and Magnetic (log).png'))
    plt.close('all')

    
def plotTstreak(anlsT, outdir='.'):
    """From anlsT, generate streak plots (field value lineout vs. time)"""
    print("Making streak plots...")
    simname = os.path.basename(os.path.normpath(outdir)) # Use the name of the output directory as a label for this sim., to use in plots etc.

    tgv = anlsT['times_ns']
    zstag = anlsT['z_stag']*1e1 # z value at stagnation, in mm, through time

    ## LINE Plots
    for lbl, rad_mm in zip(["line1", "line_500um", "line_2000um"], [0.0, 0.5, 2.0]): # Iterate over the lineouts (defined above), make plots for each
        zgv = anlsT[lbl + '_zgv'][0]*1e1 # zgv in mm
        
        fig = plt.figure(1)
        fig.clear()
        C = anlsT[lbl + '_velz'].T * 1e-5 # Velocity z component, in um/ns (converted from cm/s)
        vmax = 2*np.mean(np.abs(C))
        ax = plt.subplot(111)
        cax = ax.pcolorfast((tgv[0], tgv[-1]), (zgv[0], zgv[-1]), C, vmin=-vmax, vmax=vmax, cmap='RdBu')
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Z (mm)")
        ax.set_xlim(0, 50)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel("Velocity, Z-component (um/ns)")
        ax.set_title(r"Temporal streak of velocity-Z-component, $\rho$=" + str(rad_mm) + ' mm')
        lstag = ax.plot(tgv, zstag, lw=1, c='black', linestyle='--', label='Stagnation') # Add the stagnation z position (in mm) to plot
        ax.legend(handles=[lstag[0]], loc='upper right')
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        fig.savefig(os.path.join(outdir, "Temporal streak (" + str(rad_mm).replace('.','p') + " mm thickness) - Velocity-z.png"), dpi=300)

        fig = plt.figure(2)
        fig.clear()
        C = anlsT[lbl + '_magtot'].T * np.sqrt(4*np.pi)*1e-6 # Get magnetic field as MegaGauss (from code_magnetic)
        vmax = 1.2 * np.max(C[(zgv > 3.0) & (zgv < 5.0),:])
        ax = plt.subplot(111)
        cax = ax.pcolorfast((tgv[0], tgv[-1]), (zgv[0], zgv[-1]), C, vmin=0, vmax=vmax, cmap='viridis')
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Z (mm)")
        ax.set_xlim(0, 50)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel("Magnetic field magnitude (MGs)")
        ax.set_title(r"Temporal streak of magnetic field magnitude, $\rho$=" + str(rad_mm) + " mm")
        lstag = ax.plot(tgv, zstag, lw=1, c='black', linestyle='--', label='Stagnation') # Add the stagnation z position (in mm) to plot
        ax.legend(handles=[lstag[0]], loc='upper right')
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        fig.savefig(os.path.join(outdir, "Temporal streak (" + str(rad_mm).replace('.','p') + " mm thickness) - B magnitude.png"), dpi=300)    

        fig = plt.figure(3)
        fig.clear()
        Bsq_Tsq = (anlsT[lbl + '_magsqr'].T * 4*np.pi)*1e-8 # Magnetic field squared, in tesla squared
        EBdens_Jpmm3 = 0.5 * (Bsq_Tsq / sc.mu_0) * 1e-9
        C = EBdens_Jpmm3
        vmax = 1.2 * np.max(C[(zgv > 3.0) & (zgv < 5.0),:])
        ax = plt.subplot(111)
        cax = ax.pcolorfast((tgv[0], tgv[-1]), (zgv[0], zgv[-1]), C, vmin=0, vmax=vmax, cmap='viridis')
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Z (mm)")
        ax.set_xlim(0, 50)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel("Magnetic energy density (J/mm^3)")
        ax.set_title(r"Temporal streak of magnetic energy density, $\rho$=" + str(rad_mm) + " mm")
        lstag = ax.plot(tgv, zstag, lw=1, c='black', linestyle='--', label='Stagnation') # Add the stagnation z position (in mm) to plot
        ax.legend(handles=[lstag[0]], loc='upper right')
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        fig.savefig(os.path.join(outdir, "Temporal streak (" + str(rad_mm).replace('.','p') + " mm thickness) - Magnetic energy density.png"), dpi=300)

        fig = plt.figure(3)
        fig.clear()
        C = anlsT[lbl + '_KEdens'].T * 1e-7 * 1e-3 # Kinetic energy density, in J/mm^3
        vmax = 1.2 * np.max(C[(zgv > 3.0) & (zgv < 5.0),:])
        ax = plt.subplot(111)
        cax = ax.pcolorfast((tgv[0], tgv[-1]), (zgv[0], zgv[-1]), C, vmin=0, vmax=vmax, cmap='viridis')
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Z (mm)")
        ax.set_xlim(0, 50)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel("Kinetic energy density (J/mm^3)")
        ax.set_title(r"Temporal streak of kinetic energy density, $\rho$=" + str(rad_mm) + " mm")
        lstag = ax.plot(tgv, zstag, lw=1, c='black', linestyle='--', label='Stagnation') # Add the stagnation z position (in mm) to plot
        ax.legend(handles=[lstag[0]], loc='upper right')
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        fig.savefig(os.path.join(outdir, "Temporal streak (" + str(rad_mm).replace('.','p') + " mm thickness) - Kinetic energy density.png"), dpi=300)

        fig = plt.figure(3)
        fig.clear()
        C = anlsT[lbl + '_dens'].T # Mass density, in g/cm^3
        vmax = 1.2 * np.max(C[(zgv > 3.0) & (zgv < 5.0),:])
        ax = plt.subplot(111)
        cax = ax.pcolorfast((tgv[0], tgv[-1]), (zgv[0], zgv[-1]), C, vmin=0, vmax=vmax, cmap='viridis')
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Z (mm)")
        ax.set_xlim(0, 50)
        cbar = fig.colorbar(cax)
        cbar.ax.set_ylabel("Mass density (g/cm^3)")
        ax.set_title(r"Temporal streak of mass density, $\rho$=" + str(rad_mm) + " mm")
        lstag = ax.plot(tgv, zstag, lw=1, c='black', linestyle='--', label='Stagnation') # Add the stagnation z position (in mm) to plot
        ax.legend(handles=[lstag[0]], loc='upper right')
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
        fig.savefig(os.path.join(outdir, "Temporal streak (" + str(rad_mm).replace('.','p') + " mm thickness) - Mass density.png"), dpi=300)

    plt.close('all')
    
def plotT(anlsT, outdir='.'):
    """Custom plotting callback"""
    #Comment out what you don't want for a less complete plotting analysis
    plotTstreak(anlsT, outdir=outdir)
    plotTsph(anlsT, outdir=outdir)
    
    return anlsT

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    ## USER SIM DEFINITION
    # Note that this will run in parallel, so no prints please
    # simname # Can be anything. Will get stamped on the plots and set the name of outdir. No spaces, please.
    # datdir # Directory holding HDF5 FLASH outputs
    # basenm # Prefix name of the HDF5 FLASH outputs, ending with "_"

    simnames = [None]*3
    datdirs = [None]*3
    basenms = [None]*3
    
    simnames[1] = "OMEGA_NLUF6"
    datdirs[1] = r'/projects/Omega-NIF_Exp/tzeferac/NLUF6grpFINAL/SCRIPT5/RUN1'
    basenms[1] = r'omega2015_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file

    simnames[0] = "NIF_TDYNO_BAND"
    datdirs[0] = r'/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND'
    basenms[0] = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file

    simnames[2] = "NIF_TDYNO_200KJ"
    datdirs[2] = r'/projects/CosmicLaser/tzeferac/SCRIPT1/RUN2'
    basenms[2] = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file


    for simname, datdir, basenm in zip(simnames, datdirs, basenms):
        ## READ IN THE SIM NAMES
        fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_????') # yt-time-series filename pattern for plot files
        outroot=r'/home/sfeister/myouts/'
        
        if rank == 0:
            print("STARTING WORK ON: " + simname)
            #outdir = sf.subdir(outroot, "MagAnalysis-{:%Y-%m-%d_%H%M}".format(datetime.now()))
            outdir = sf.subdir(outroot, simname)
            #outdir = sf.subdir(outroot,"scratch_dev")
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

