#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
magtrack3.py: Analyze magnetic fields, tracking the zero-velocity front

Created by Scott Feister on Fri Oct 19 2016

USAGE:
mpirun -np 12 python magtrack3.py

CHANGELOG:
2016-10-19 Spun off from "tstest5.py", which was to test parallel processing
2016-11-18 Changed from magtrack2.py to magtrack3.py
           Removed "stagnation" label
           Added more variables
2016-11-29 Added a bunch of plasma variables in morefields.py, then included here.
           Got rid of slow units determination
2016-06-28 Spun off from magtrack3.py; for use in cshoc tests
           
TODO:
* Organize the calls for lineout analysis of various widths
* Clean up fatline function (from flyt) to reduce memory overhead of multiple calls
* Incorporate a kinetic energy vs. magnetic energy analysis
* Better incorporate simulation name? Not sure.
* Clean up the lineouts plotting outputs
* Update magnetic field energy so it's not vacuum energy (mu_0) but the real mu?
* Sort out the use of included derived fields where appropriate
"""


# SHOULD BE FALSE unless this is the most simple re-analysis
quick = False # If true, assume analysis has already been completed and only remake plots

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
    
    anlsD['myflds'] = ['nion', 'nele', 'coulog', 'Pm', 'RmperL', 'ReperL', 'tele', 'tion', 'iimfp', 'abar', 'zbar'] # A more general list of fields
    fldids = anlsD['myflds'] + ['velz', 'magtot', 'dens', 'KEdens', 'magsqr', 'cs', 'mach']

    #fldids = ['velz', 'mach'] #TEST
    
    # Extract a lineout at X=Y=0 for velocity
    anlsD['line1_zgv'], fldgv = flyt.lineout2(ds, flds=fldids, axis=2, npts=2000)
    for fldid in fldids:
        anlsD['line1_' + fldid] = fldgv[fldid]
        anlsD['line1_' + fldid + '_units'] = str(fldgv[fldid].units)

    ## TODO: Organize this fatter and fatter lineout section. Also, fix the fatline function so it takes all the fields at once (minimize overhead)
    # Extract a fatter lineout at X=Y=0 for similar quantities
    rad_mm = 0.5 # Radius of the lineout cylinder
    anlsD['line_500um_zgv'], fldgv = flyt.fatline2(ds, rad_mm, flds=fldids, npts=2000)
    for fldid in fldids:
        anlsD['line_500um_' + fldid] = fldgv[fldid]
        anlsD['line_500um_' + fldid + '_units'] = str(fldgv[fldid].units)
        
    # Extract an even fatter lineout at X=Y=0 for similar quantities
    rad_mm = 1.0 # Radius of the lineout cylinder
    anlsD['line_1000um_zgv'], fldgv = flyt.fatline2(ds, rad_mm, flds=fldids, npts=2000)
    for fldid in fldids:
        anlsD['line_1000um_' + fldid] = fldgv[fldid]
        anlsD['line_1000um_' + fldid + '_units'] = str(fldgv[fldid].units)

    # Extract a yet-even fatter lineout at X=Y=0 for similar quantities
    rad_mm = 2.0 # Radius of the lineout cylinder
    anlsD['line_2000um_zgv'], fldgv = flyt.fatline2(ds, rad_mm, flds=fldids, npts=2000)
    for fldid in fldids:
        anlsD['line_2000um_' + fldid] = fldgv[fldid]
        anlsD['line_2000um_' + fldid + '_units'] = str(fldgv[fldid].units)

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

def customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, cmap=None, vmin=None, vmax=None, logit=False):
    """ One-off, custom plotting function for Scott's use in streak plots 
    Creates and saves a figure based on specific and inflexible inputs """
    fig = plt.figure(7)
    fig.clear()
    if logit: # Logarithm of (data * fld_mult)
        C = np.log10(anlsT[lbl + '_' + fld].T * fld_mult) # Read in the data, along with a multiplier
    else: # Normal; (data * fld_mult)
        C = anlsT[lbl + '_' + fld].T * fld_mult # Read in the data, along with a multiplier
        
    if vmax is None:
        #vmax = 1.2 * np.max(np.abs(C[(zgv > 3.0) & (zgv < 5.0),:])) # TDYNO-style
        vmax = 1.2 * np.max(np.abs(C[(zgv > -2.0) & (zgv < 2.0),:])) # CSHOC-style
    
    if vmin is None:
        if np.min(C) < 0: # Do we plot symmetric about zero, or all greater than zero?
            vmin = -vmax
            if cmap is None:
                cmap = 'RdBu' # Diverging colormap
        else:
            vmin = 0
    if cmap is None:
        cmap = 'viridis' # Sequential colormap
        
    if units is None: # Get the units if not specified (units = None)
        units = anlsT[lbl + '_' + fld + '_units'][0]
    ax = plt.subplot(111)
    cax = ax.pcolorfast((tgv[0], tgv[-1]), (zgv[0], zgv[-1]), C, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Z (mm)")
    ax.set_xlim(0, 10) # TEMPORARY TEMPORAL LIMITS
    #ax.set_ylim(2, 6) # TEMPORARY Z LIMIT ENFORCEMENT
    cbar = fig.colorbar(cax)
    if units == "dimensionless":
        ylabel = name
    else:
        ylabel = name + " (" + units + ")"
    cbar.ax.set_ylabel(ylabel)
    ax.set_title(r"Temporal streak of " + name + r", $\rho$=" + str(rad_mm) + " mm")
    fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower-right footer
    fig.savefig(os.path.join(outdir, "Temporal streak (" + str(rad_mm).replace('.','p') + " mm thickness) - " + name + ".png"), dpi=300)
    return fig
        
def plotTstreak(anlsT, outdir='.'):
    """From anlsT, generate streak plots (field value lineout vs. time)"""
    print("Making streak plots...")
    simname = os.path.basename(os.path.normpath(outdir)) # Use the name of the output directory as a label for this sim., to use in plots etc.

    tgv = anlsT['times_ns']
    zstag = anlsT['z_stag']*1e1 # z value at stagnation, in mm, through time

    ## LINE Plots
    for lbl, rad_mm in zip(["line1", "line_500um", "line_1000um", "line_2000um"], [0.0, 0.5, 1.0, 2.0]): # Iterate over the lineouts (defined above), make plots for each
        zgv = anlsT[lbl + '_zgv'][0]*1e1 # zgv in mm
        
        fld = 'mach' # Mach number (local fluid velocity / local sound speed)
        fld_mult = 1 # Amount to multiply by
        name = "Mach number"
        units = "dimensionless"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'velz' # Velocity z component
        fld_mult = 1e-5 # convert to um/ns (converted from cm/s) 
        name = "Velocity, Z-component"
        units = "um/ns"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'magtot' # Magnetic field magnitude
        fld_mult = np.sqrt(4*np.pi)*1e-6 # Converted into MegaGauss (from code_magnetic)
        name = "Magnetic field magnitude" # Or, "B field magnitude"
        units = "MGs"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'magsqr' # Magnetic field magnitude, squared
        fld_mult = (4 * np.pi * 1e-8) * (0.5 * 1e-9 / sc.mu_0)  # Converted first into tesla^2, then into magnetic energy density, J/mm3 (from code_magnetic**2)
        name = "Magnetic energy density" # Or, "B field magnitude"
        units = "J/mm^3"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'KEdens' # Kinetic energy density
        fld_mult = 1
        name = "Kinetic energy density"
        units = "J/mm^3"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'dens' # Mass density, in g/cm^3
        fld_mult = 1
        name = "Mass density"
        units = "g/cm^3"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'cs' # Sound speed
        fld_mult = 1
        name = "Sound speed"
        units = None # Use default units
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'nion' # Ion density
        fld_mult = 1
        name = "Ion density"
        units = "1/cm^3" # Use default units
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'nele' # Electron density (per cm^3)
        fld_mult = 1
        name = "Electron density"
        units = "1/cm^3" # Use default units
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'tion' # Ion temperature (in K)
        fld_mult = 1/11604.525 # Convert to eV
        name = "Ion temperature"
        units = "eV"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'tele' # Electron temperature (in K)
        fld_mult = 1/11604.525 # Convert to eV
        name = "Electron temperature"
        units = "eV"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'iimfp' # Ion-ion MFP
        fld_mult = 1e4 # Convert to microns
        name = "Ion-ion mean free path"
        units = "um"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)

        fld = 'RmperL' # Magnetic Reynold's number, per L (in cm)
        fld_mult = 0.06 # Convert to Magnetic Reynold's number, assume L = 600 um (0.06 cm)
        name = "Magnetic Reynolds Number"
        units = "dimensionless"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, vmin=0)

        fld = 'abar' # Average atomic weight (amu)
        fld_mult = 1
        name = "Average atomic weight"
        units = "a.m.u."
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, cmap='Paired')

        fld = 'zbar' # Average ionization
        fld_mult = 1
        name = "Average ionization"
        units = "dimensionless"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, cmap='Paired')

        fld = 'Pm' # Prandtl number
        fld_mult = 1
        name = "Log10[Prandtl number]"
        units = "dimensionless" # Use default units
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, vmin=-5, vmax=7, cmap='Paired', logit=True)

        fld = 'ReperL' # Reynold's number, per L (in cm)
        fld_mult = 0.06 # Convert to Reynold's number, assume L = 600 um (0.06 cm)
        name = "Reynolds Number"
        units = "dimensionless"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, vmin=0)

        fld = 'coulog' # Coulomb logarithm
        fld_mult = 1
        name = "Coulomb logarithm"
        units = "dimensionless"
        fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname, vmin=0)

        #for fldid in anlsT['myflds'][0]: # Iterate over all the junky fields
        #    fld = fldid
        #    fld_mult = 1
        #    name = fldid
        #    units = None # Use default units
        #    fig = customstreak1(anlsT, lbl, fld, fld_mult, name, units, zgv, tgv, rad_mm, outdir, simname)       
        
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

    nsims = 1
    
    simnames = [None]*nsims
    datdirs = [None]*nsims
    basenms = [None]*nsims

    simnames[0] = "CSHOC"
    datdirs[0] = r'/project/tzeferacos/cshoc_outs/run-demo2_FL4prod_mres_20170620'
    basenms[0] = r'cshoc_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
    
    #simnames[2] = "OMEGA_NLUF4"
    #datdirs[2] = r'/projects/Omega-NIF_Exp/tzeferac/NLUF6grpFINAL/SCRIPT4/RUN1'
    #basenms[2] = r'omega2015_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
    
    #simnames[3] = "OMEGA_NLUF6"
    #datdirs[3] = r'/projects/Omega-NIF_Exp/tzeferac/NLUF6grpFINAL/SCRIPT5/RUN1' # BROKEN LINK?
    #basenms[3] = r'omega2015_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file

    #simnames[1] = "NIF_TDYNO_BAND"
    #datdirs[1] = r'/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND'
    #basenms[1] = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file

    #simnames[0] = "NIF_TDYNO_200KJ"
    #datdirs[0] = r'/projects/CosmicLaser/tzeferac/SCRIPT1/RUN3'
    #basenms[0] = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file


    for simname, datdir, basenm in zip(simnames, datdirs, basenms):
        ## READ IN THE SIM NAMES
        fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_????') # yt-time-series filename pattern for plot files
        # hdf5_plt_cnt_???[0,2,4,6,8]
        # hdf5_plt_cnt_???[0,5]
        # hdf5_plt_cnt_??[0,2,4,6,8]0
        # 'hdf5_plt_cnt_??[0,5]0'
        #outroot=r'/home/sfeister/myouts/'
        #outroot=r'/home/sfeister/myouts/scratch_dev' # Cooley
        outroot=r'/project/tzeferacos/cshoc_outs/Analyzed' # Midway2
        
        
        if rank == 0:
            print("STARTING WORK ON: " + simname)
            #outdir = sf.subdir(outroot, "MagAnalysis-{:%Y-%m-%d_%H%M}".format(datetime.now()))
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

