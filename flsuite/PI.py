#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
PI.py: Proton imaging tools; read and make plots of proton imaging outputs made by FLASH

Created by Scott Feister on Fri Oct 07 13:36:27 2016

This module is rather fragile, as it depends on the exact output syntax of the Proton Imaging module as of FLASH 4.3

CHANGELOG:
2016-10-07 Created piIO.py to focus mainly on I/O from proton imaging
2016-10-07 Renamed to PI.py, to account for lots of proton plotting tools absorbed from piPlot.py

TODO:
* Implement checks on the file format syntaxes to avoid plowing forward if file formats are wrong.
"""

import os
# Stuff for plotting parts
import matplotlib
matplotlib.use('Agg') # Display choice "Agg" for headless servers

import re
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches # Using this example to draw circles: http://matplotlib.org/examples/shapes_and_collections/artist_reference.html
import sftools as sf


# Works only for single-beam radiography!! TODO: Implement error if has more than one beam
def piRead(fn):
    """Read ProtonImagingMainPrint, ProtonDetectorsPrint, or ProtonBeamsPrint as key, value pairs into a dict file"""
    with open(fn, 'r') as f:
        buff = f.read()

    d = {}
    p = re.compile('^\s*(.*?)\s*=\s*(.*?)\s*$', re.MULTILINE)
    matches = p.findall(buff)
    for m in matches: # TODO: Split into appropriate sections!
        d[m[0]] = m[1]
        
    return d

# Expect this to be replaced with something more elegant soon.
def piHugeAnalysis(PIdir, basenm=r"tdyno2016PI_", simname=None, outdir=None, pitdiam_um = 10, bin_um = 332.6):
    """Performs a massive, custom analysis. Outputs plots in PIdir
    
    PIdir: Path to folder containing the PI outputs like blahblah_ProtonImagingMainPrint    
    simname: (String) Name user gives this sim for plotting purposes, can be anything here. If left as none, defaults to basename.
    pitdiam_um = 10 # Pit diameter, in microns
    bin_um = 332.6 # Length of one side of the square CR39 bin, in microns
    
    TODO: Deprecate this function. Tries to do too much, not modular, too disorganized.
    """

    if not simname:
        simname=basenm
    
    if not outdir:
        outdir=PIdir
    
    ################## GENERAL ANALYSIS ####################
     # CR39 pit and readout properties
     # Pit diameter, in microns

    ## All of this could somehow be packaged into something higher-level? Like, maybe a reader function? Not sure...
    # Read in a few files    
    #maindict = piRead(os.path.join(PIdir, basenm + r"ProtonImagingMainPrint.txt"))
    try:
        detdict = piRead(os.path.join(PIdir, basenm + r"ProtonDetectorsPrint.txt"))
    except:
        detdict = piRead(os.path.join(PIdir, basenm + r"ProtonImagingDetectors.txt"))
    beamdict = piRead(os.path.join(PIdir, basenm + r"ProtonBeamsPrint.txt"))

    # Get some important variables from the read-ins
    dist_cm = float(detdict["Detector distance from beam capsule center"]) # Distance from capsule center to CR39 center
    width_cm = float(detdict["Detector square side length (cm)"]) # Length of one side of the square CR39 detector, in cm
    apdegs = np.around(np.rad2deg(float(beamdict["Beam aperture angle (rad)"])),11) # Cone aperture angle (full-angle), in degrees. Round to 11th decimal place (to avoid rad2deg mismatches; present 1.0 rather than 0.99999999999980005)
    protMeV = float(beamdict["Proton energy (in MeV)"]) # Proton energy, in MeV
    
    ## Extract timestamp from filenames
    fns = sf.getfns(PIdir, prefix = basenm + 'ProtonDetectorFile')

    # Loop over all the functions
    for fn in fns:
        p = re.compile(basenm + r'ProtonDetectorFile([0-9]+)_(\S*)') # Strip timestamp off filename end, e.g. tdyno2016PI_ProtonDetectorFile01_2.200E-08 ==> 2.2000E-08
        m = p.findall(fn)
        #detnum = int(m[0][0]) # Detector ID number (e.g. 1, 2, 3,..)
        time_ns = float(m[0][1])*1e9 # Time step in nanoseconds
        tlabel = str(m[0][1])
        
        # Read in the proton file
        dat = np.genfromtxt(fn)
        if len(np.atleast_1d(dat.flatten())) < 1:
            print("File contents empty : " + fn + ". Moving on...")
            continue
        
        xy_cm = (dat[:,(0,1)] - 0.5) * width_cm # Convert scatter points from 0 to 1 grid up to centimeters
    
        ##################### DETAILED ANALYSIS ###################
        ## Make some plots
        # Calculate solid angle of CR39
        alph = width_cm / (2*dist_cm)
        sr = 4 * np.arccos(np.sqrt( (1 + 2 * alph**2) / (1 + alph**2)**2 ) ) # CR39 solid angle relative to capsule center, in s.r.
        
        # Calculate the undeflected beam radius
        protrad_cm = dist_cm * np.tan(np.deg2rad(apdegs/2))# Radius of undeflected cone on CR39, in centimeters
    
        print "Histogramming..."
        # Bin the data, according to the square bin edge size
        # Following example at: http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram2d.html
        bins_cm = np.arange(-width_cm/2, width_cm/2, bin_um*1e-4) # 1D array of bin edges, in centimetres
        H, xedges, yedges = np.histogram2d(xy_cm[:,0], xy_cm[:,1], bins=bins_cm)
        
        print "Making PI plot..."
        ## Figure 1: Main radiograph
        fig = plt.figure(1)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_title('$FLASH\ protons:$ ' + simname + ', ' + str(apdegs) + '$^\circ$ ap.')
        X, Y = np.meshgrid(xedges, yedges)
        vmax = np.ceil(np.max(H)/5.0)*5.0 # Round up to nearest 5 for the colormap max
        #vmax = 150.0
        cax = ax.pcolormesh(X, Y, H.T, cmap='Greys', vmin=0, vmax=vmax) # Transpose needed because H array is organized H[xindex, yindex] but this is flipped from what pcolormesh, meshgrid output. (E.g. X[:,1] gives a uniform number)
        # Draw a circle, also
        circle = mpatches.Circle((0,0), radius=protrad_cm, fill=False, edgecolor="blue", linestyle="--", label='Undeflected')
        ax.add_patch(circle)
        ax.set_xlim([np.min(bins_cm), np.max(bins_cm)])
        ax.set_ylim([np.min(bins_cm), np.max(bins_cm)])
        ax.set_aspect('equal')
        plt.legend()
    
        # Add colorbar, make sure to specify tick locations to match desired ticklabels
        plt.xlabel('CR39, X (cm)')
        plt.ylabel('CR39, Y (cm)')
        #plt.colorbar(label='Protons/bin')
        cbar = fig.colorbar(cax, label='Protons/bin')
        plt.tight_layout()
        #sb.jointplot(dat[:,0], dat[:,1], kind='hex')
    
        tstring =  't=' + "{:.1f}".format(time_ns) + " ns" # Time string
        Estring =  "{:.1f}".format(protMeV) + " MeV" # Proton energy string
        ax.text(0.05, 0.95, tstring, fontsize=18, color='black', transform=ax.transAxes, horizontalalignment='left', verticalalignment='top') # Upper left within axis (transform=ax.transAxes sets it into axis units 0 to 1)
        ax.text(0.05, 0.03, Estring, fontsize=24, color='maroon', transform=ax.transAxes, horizontalalignment='left', verticalalignment='bottom') # Lower left within axis
        fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower right in figure units
        plt.savefig(os.path.join(outdir, "Radiograph_" + tlabel + ".png"), dpi=300)
        
        
        if dat.shape[1] >= 5:
            print("Making Energy spectrum plot")
            velxyz_mps = dat[:,(3,4,5)] * 1e-2 # Fourth, fifth, sixth columns are vx, vy, vz in cm/s. Convert to m/s
            vel_mps = np.sqrt(np.sum(velxyz_mps**2, 1)) # Magnitude of velocity, in m/s
            print("Shape of velxyz_mps: ", velxyz_mps.shape)
            print("Shape of vel_mps: ", vel_mps.shape)
            print("Max of velxyz_mps:", velxyz_mps.max())
            print("Min of velxyz_mps:", velxyz_mps.min())
            print("NaN in velxyz_mps:", np.isnan(velxyz_mps).sum())

            gamma = 1/np.sqrt(1 - (vel_mps**2 / sc.c**2)) # Relativistic gamma factor
            print("Mean gamma factor: " + str(np.mean(gamma)))
            KE_J = (gamma - 1) * sc.m_p * sc.c**2  # Relativistic energy, assuming proton mass
            KE_MeV = (KE_J / sc.e) * 1e-6
            print("Shape of KE_MeV: ", KE_MeV.shape)

            histrange = protMeV + np.array([-0.25, 0.25]) # Range for the histogram, in MeV
            fig = plt.figure(2)
            plt.clf()
            n, bins, patches = plt.hist(KE_MeV, 100, range=histrange, normed=True, facecolor='green', alpha=0.75)
            plt.title("Final proton spectrum")
            plt.xlabel('Energy (MeV)')
            plt.ylabel('Number (a.u.)')
            fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower right in figure units
            plt.savefig(os.path.join(outdir, "ESpec_" + tlabel + ".png"), dpi=300)

            fig = plt.figure(3)
            plt.clf()
            n, bins, patches = plt.hist(KE_MeV, 100, normed=True, facecolor='green', alpha=0.75)
            plt.title("Final proton spectrum (auto-ranged)")
            plt.xlabel('Energy (MeV)')
            plt.ylabel('Number (a.u.)')
            fig.text(0.99, 0.01, simname, horizontalalignment='right') # Lower right in figure units
            plt.savefig(os.path.join(outdir, "ESpecAuto_" + tlabel + ".png"), dpi=300)
        else:
            print("No velocity data found; perhaps this was a pre-Oct2016 version of FLASH proton imaging.")
#        ## Figure 2 & 3: Other stuff
#        #TODO: Plot the densest cell in CR-39 fashion??
#        
#        # Get the index of the densest histogrammed cell
#        [imax, jmax] = np.unravel_index(H.argmax(), H.shape)
#        xmin = xedges[imax] # Edge of the bin
#        xmax = xedges[imax + 1]
#        ymin = yedges[jmax]
#        ymax = yedges[jmax + 1]
#        
#        ct = (dat_cm[:,0] > xmin) & (dat_cm[:,0] < xmax) & (dat_cm[:,1] > ymin) & (dat_cm[:,1] < ymax)
#        subdat_cm = dat_cm[ct,:]
#        numpits = subdat_cm.shape[0] # Total number of pits to be plotted
#        print "Making the pits plot (" + str(numpits) + " pits)..."
#        
#        
#        fig = plt.figure(2)
#        ax = plt.subplot(111)
#        pitrad_cm = (pitdiam_um / 2) * 1e-4 # Pit diameter converted to centimeters
#        for i in range(numpits):
#        	spot = mpatches.Circle((subdat_cm[i,0], subdat_cm[i,1]), radius=pitrad_cm, fill=True, edgecolor="blue", linestyle="-")
#        	ax.add_patch(spot)
#        
#        plt.title("CR39 bin with the most protons")
#        plt.xlim(xmin, xmax)
#        plt.ylim(ymin, ymax)
#        plt.xlabel('CR39, X (cm)')
#        plt.ylabel('CR39, Y (cm)')
#        
#        #colors = 0.5 * np.ones(subdat_cm.shape[0])
#        #area = 5
#        #plt.scatter(subdat_cm[:,0], subdat_cm[:,1], s=area, c=colors, alpha=0.5)
#        #plt.scatter([0], [0], s=area, c=colors, alpha=0.5)
#        
#        ax.set_aspect('equal')
#        plt.savefig(os.path.join(PIdir, "SamplePits.png"))
        print("Done with time: " + tlabel)
        
        # Big ole scatter plot of all pits
        #fig = plt.figure(3)
        #plt.scatter(dat_cm[:,0], dat_cm[:,1], s=1, alpha=0.1)
        #plt.xlabel("X")
        #plt.ylabel("Y")
        #plt.savefig(os.path.join(PIdir, "AllPits.png"))
    
    return 0
    
if __name__ == "__main__":
    ## USER PARAMETERS
    PIdirs = [x[0] for x in os.walk(r"C:\Users\Scott\Documents\temp\TDYNO vids\PI junk\PINHOLE")][1:]
    #PIdir = r"C:\Users\Scott\Documents\temp\TDYNO vids\PI junk\PINHOLE\PIN-10-1p0-15"
    for PIdir in PIdirs:
        print PIdir
        if not os.path.isfile(os.path.join(PIdir, "Radiograph.png")):
            piHugeAnalysis(PIdir, basenm=r"tdyno2016_", simname='TDYNO_BAND')
        else:
            print "Skipping..."
