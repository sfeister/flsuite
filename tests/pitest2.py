#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
pitest2.py: Test of reading a proton imaging plotfile contents

Created by Scott Feister on Wed Dec 07 15:31:48 2016

Makes some very nice 3D plots.

Changelog:
2016-12-12 Spun off this program from trimtest.py (analysis of TRIM data)

TODO:
* Better framework for selection of "good" particles
* Pre-calculate a variety of values? Keep velocity, position, start/end, as vectors?

"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

if __name__ == "__main__":
    ## Load into working memory and sort
    fn = r"C:\Users\Scott\Documents\temp\subdata.hdf5"
    f = h5py.File(fn, 'r')
    dat = f["ProtonData"][...]
    ix_arr = np.argsort(dat[:,0], kind='mergesort') # Sort by particle ID # (index 0). Mergesort keeps the paths in order, for a given ID # # TODO: Sort in batch sizes?
    dat = dat[ix_arr, :] # Sorted array
    
    # Split into sub-arrays / tracks
    # http://stackoverflow.com/questions/19125661/find-index-where-elements-change-value-numpy
    idarr = dat[:, 0] # Array of particle ID numbers
    splitix = np.where(idarr[:-1] != idarr[1:])[0] + 1 # Array of indices at which particle ID changes (split here)
    edgeix = np.hstack((0, splitix, len(idarr))) # Array of start and stop indices; N + 1 long, where N is number of particle tracks
    traklens = np.diff(edgeix) # Array of particle track lengths
    ntraks = len(traklens) # Number of particle tracks
    
    ## Option 1: Split into a list of arrays (Not used, for now)
    tracks = np.split(dat, splitix, axis=0) # Original array, split up into discrete lists by particle
    #plt.figure(1)
    #plt.clf()
    
    ## Option 2: Split into a larger array with one more dimension
    tlmax = 601 # Max track length (number of steps); Either the max number in the arrays, or 50 steps
    trakarr = np.empty((ntraks, tlmax, dat.shape[1]))
    trakarr.fill(np.nan)
    
    for i in range(ntraks):
        trakarr[i,:min(traklens[i], tlmax),:] = dat[edgeix[i]:edgeix[i] + min(traklens[i], tlmax),:]
    
    dtrakarr = np.diff(trakarr, axis=1) # Difference along the traveling axis
    dr = np.sqrt(dtrakarr[:,:,1]**2 + dtrakarr[:,:,2]**2 + dtrakarr[:,:,3]**2) # N x Steps-1 long, Path length change between steps
    pathl = np.concatenate((np.zeros((dr.shape[0],1)), np.cumsum(dr, axis=1)), axis=1) # N x Steps long, Total path length up to this step

    posxyz = np.rollaxis(trakarr[:,:,1:],2,0) # Positions along the trajectory (Component x Particle x Step#)
    velxyz = np.gradient(posxyz, axis=2) # Non-normalized velocities along the trajectories (Component x Particle x Step#)
    velxyz = velxyz / np.sqrt(np.sum(velxyz**2, axis=0)) # Normalized velocities along the trajectories, normalized to 1 (Component x Particle x Step#)

    velmag = np.sqrt(np.sum(velxyz**2, axis=0)) # This is actually always 1, in our case
    theta = np.arccos(velxyz[1]/velmag)
    
    accxyz = np.gradient(velxyz, axis=2) # Accelerations, using a normalized velocity
    accmag = np.sqrt(np.sum(accxyz**2, axis=0)) # Magnitude of acceleration, using the normalized velocity
    
    fig = plt.figure(1)
    fig.clear()
    ax = fig.add_subplot(111)

    ntplot = 1000 # Number of tracks to plot
    maxpathl = 100.0
    
    ct = pathl[:ntplot,:] > maxpathl # nan condition
    xvals = trakarr[:ntplot,:,1] # X values along trajectory
    yvals = trakarr[:ntplot,:,2]
    zvals = trakarr[:ntplot,:,3]
    lvals = pathl[:ntplot,:] # Path length values along trajectory
    xvals[ct] = np.nan
    yvals[ct] = np.nan
    zvals[ct] = np.nan
    lvals[ct] = np.nan

    # (JANKY) Extract X,Y,Z endpoints
    xend = np.zeros((xvals.shape[0])) # 1D array to hold endpoints in X
    yend = np.zeros_like(xend)
    zend = np.zeros_like(xend)
    tof = np.zeros_like(xend) # 1D array of time of flights
    ixend = np.nanargmax(yvals,1)
    for i in range(len(ixend)):
        j = ixend[i]
        xend[i] = xvals[i, j]
        yend[i] = yvals[i, j]
        zend[i] = zvals[i, j]
        tof[i] = lvals[i, j] # TODO: Sort by tof?
    
    scl = velxyz[1,:,0] / np.sqrt(np.sum(velxyz[:,:,0]**2, axis=0)) # Adjustment for TOF; Y-component of init. velocity
    tofadj = tof*scl # Time of flight, adjusted for starting velocity angles

    ## 2D head-on plot of particle trajectories
    fig = plt.figure(1)
    fig.clear()
    ax = fig.add_subplot(111)
    ax.plot(xvals.T, zvals.T)
    ax.set_aspect('equal')
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    ax.set_title("Head-on view of particle tracks")
    #ax.set_ylim(-3000, 3000)
    #ax.set_xlim(0.27, 0.29)
    
    ## 3D view of particle trajectories
    fig = plt.figure(2)
    fig.clear()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(xvals.shape[0]):
        ax.plot(xvals[i,:], yvals[i,:], zvals[i,:], alpha=0.1)
    ax.scatter(xvals[:,0], yvals[:,0], zvals[:,0], marker='.', c='k', alpha=0.15, s=10) # Initial posits
    ax.scatter(xend, yend, zend, marker='.', c='k', alpha=0.15, s=10) # Initial posits
    ax.set_aspect('equal')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.elev=10
    ax.azim=40
    ax.set_title("3D view of particle tracks")

    ## Histogram of times of flight, maps of init and final posits
    fig = plt.figure(3)
    fig.clear()

    ax0 = fig.add_subplot(221)
    ax0.hist(tof, bins=100,range=(0.62, 0.64))
    ax0.set_xlim(0.62, 0.64)
    ax0.set_xlabel("TOF (a.u.)")
    ax0.set_ylabel("Number density (a.u.)")
    ax0.set_title("Time of flight")
    ax0.set_yscale('log')

    ax1 = fig.add_subplot(222)
    ax1.hist(tofadj, bins=100,range=(0.62, 0.64))
    ax1.set_xlim(0.62, 0.64)
    ax1.set_xlabel(r"TOF (adjusted for $\theta_i$, a.u.)")
    ax1.set_ylabel("Number density (a.u.)")
    ax1.set_title("Adjusted time of flight")
    ax1.set_yscale('log')
    
    ax2 = fig.add_subplot(223)
    ax2.hexbin(xvals[:,0], zvals[:,0], gridsize=20, extent=(-0.3,0.3,0.1,0.7))
    ax2.set_xlabel("X (init.)")
    ax2.set_ylabel("Z (init.)")
    ax2.set_aspect('equal')
    ax2.set_title("Initial positions")

    ax3 = fig.add_subplot(224)
    ax3.hexbin(xend, zend, gridsize=20, extent=(-0.3,0.3,0.1,0.7))
    ax3.set_xlabel("X (final)")
    ax3.set_ylabel("Z (final)")
    ax3.set_aspect('equal')
    ax3.set_title("Final positions")

    plt.tight_layout()
    
    ## 3D view of only the best particle trajectories
    ix_tofsrt = np.argsort(tofadj)[::-1] # A sorting by time-of-flight
    slct = np.sort(ix_tofsrt[:25]) # Selected indices for plotting, in ascending order
    
    fig = plt.figure(4)
    fig.clear()
    ax = fig.add_subplot(111)
    ax.plot(xvals[slct,:].T, zvals[slct,:].T)
    ax.set_aspect('equal')
    ax.set_xlim(-0.3, 0.3)
    ax.set_ylim(0.1, 0.7)
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    ax.set_title("Head-on view of TOP-25 particle tracks")

    fig = plt.figure(5)
    fig.clear()
    ax = fig.add_subplot(111, projection='3d')
    for i in slct:
        ax.plot(xvals[i,:], yvals[i,:], zvals[i,:], alpha=0.1)
    ax.scatter(xvals[:,0], yvals[:,0], zvals[:,0], marker='.', c='k', alpha=0.05, s=10) # Initial posits
    ax.scatter(xend, yend, zend, marker='.', c='k', alpha=0.05, s=10) # Final posits
    ax.scatter(xvals[slct,0], yvals[slct,0], zvals[slct,0], marker='.', c='k', alpha=0.15, s=10) # Select Initial posits
    ax.scatter(xend[slct], yend[slct], zend[slct], marker='.', c='k', alpha=0.15, s=10) # Select Final posits
    ax.set_aspect('equal')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.elev=4
    ax.azim=43
    ax.set_title("3D view of TOP-25 particle tracks")
    
    fig = plt.figure(8)
    fig.clear()
    ax1 = fig.add_subplot(111)
    ax1.plot(np.nanmean(accmag, 0))
    ax1.set_title("Acceleration vs. Particle step (All particles)")
    ax1.set_ylabel("Acceleration, averaged over particles (a.u.)")
    ax1.set_xlabel("Particle step")

    fig = plt.figure(9)
    fig.clear()
    ax2 = fig.add_subplot(111)
    offset = np.arange(accmag.shape[0])*0.004
    ax2.plot(posxyz[1,:100,:].T, (accmag[:100,:].T - offset[:100]), lw=0.5)
    ax2.set_title("Acceleration vs. Y depth (100 particles)", fontsize=20)
    ax2.set_xlabel("Spatial Y", fontsize=20)
    ax2.set_ylabel("Acceleration + offset (a.u.)", fontsize=20)
    
    
    fig = plt.figure(10)
    fig.clear()
    ax = fig.add_subplot(221)
    theta_rel = (theta.T - theta[:,0].T).T
    dtheta = np.gradient(theta, axis=1)
    ax.pcolormesh(theta_rel[:,:200])
    
    ax2 = fig.add_subplot(222)
    ax2.plot(posxyz[1,:10,:].T, np.rad2deg(theta[:10,:].T))
    ax2.set_title("Theta vs. Y depth (10 particles)")
    ax2.set_ylabel("Theta (degs)")
    ax2.set_xlabel("Spatial Y")

    ax3 = fig.add_subplot(223)
    ax3.plot(dtheta[:2,:].T)
    ax3.set_ylim(-0.006, 0.006)
    ax3.set_title("dTheta vs. Step # (2 particles)")
    ax3.set_ylabel("dTheta (drad)")

    ax4 = fig.add_subplot(224)
    ax4.plot(pathl[:10,:].T, np.rad2deg(theta[:10,:].T))
    ax4.set_title("Theta vs. Path length (10 particles)")
    ax4.set_ylabel("Theta (degs)")
    ax4.set_xlabel("Path length")
    

    ## Range of x particle plot movie (uncomment savefig line)
    xc = 0.05
    zc = 0.45
    fig = plt.figure(11)
    fig.clear()
    gs = gridspec.GridSpec(2, 2)

    ax0 = fig.add_subplot(gs[0,1])
    ax1 = fig.add_subplot(gs[1,1])
    ax2 = fig.add_subplot(gs[:,0], projection='3d')

    i = 0
    for xc in np.linspace(-0.15, 0.15, 10):
        for zc in (0.4 + np.linspace(-0.15, 0.15, 10)):
            ct2 = np.sqrt((xc - xvals[:,0])**2 + (zc - zvals[:,0])**2) < 0.03
            
            
            ax0.clear()
            ax0.scatter(xvals[ct2,0], zvals[ct2,0], marker='.', alpha=0.5, s=10)
            ax0.scatter(xvals[:,0], zvals[:,0], marker='.', alpha=0.1, s=10)
            ax0.set_title("Initial positions")
            ax0.set_aspect('equal')
            
            ax1.clear()
            ax1.scatter(xend[ct2], zend[ct2], marker='.', alpha=0.5, s=10)
            ax1.scatter(xend, zend, marker='.', alpha=0.1, s=10)
            ax1.set_title("Final positions")
            ax1.set_aspect('equal')
            
            ax2.clear()
            ax2.scatter(xvals[:,0], yvals[:,0], zvals[:,0], marker='.', c='k', alpha=0.05, s=10) # Initial posits
            ax2.scatter(xend, yend, zend, marker='.', c='k', alpha=0.05, s=10) # Final posits
            
            for j in np.where(ct2)[0]:
                ax2.plot(xvals[j,:], yvals[j,:], zvals[j,:], alpha=0.2)
            ax2.scatter(xvals[ct2,0], yvals[ct2,0], zvals[ct2,0], marker='.', c='k', alpha=0.15, s=10) # Select Initial posits
            ax2.scatter(xend[ct2], yend[ct2], zend[ct2], marker='.', c='k', alpha=0.15, s=10) # Select Final posits
            ax2.set_aspect('equal')
            ax2.set_xlabel("X")
            ax2.set_ylabel("Y")
            ax2.set_zlabel("Z")
            ax2.elev=4
            ax2.azim=33
            ax2.set_title("3D view of " + str(np.sum(ct2)) + " particles")
            
            if i == 0:
                plt.tight_layout()

            #fig.savefig(os.path.join(r'C:\Users\Scott\Documents\TDYNO\plots\2016-12-16 raster movie 1', str(i).zfill(5) + '.png'))
            i+=1
    
    ## Individual particle plot movie
    fig = plt.figure(12)
    fig.clear()
    ax = fig.add_subplot(121, projection='3d')
    #ax2 = fig.add_subplot(212)
    ax3 = fig.add_subplot(122)
    #ax4 = fig.add_subplot(224)

#    for i in ix_tofsrt:
    for i in []:
        ax.clear()
        ax.plot(xvals[i,:], yvals[i,:], zvals[i,:], alpha=0.5)
        ax.scatter(xvals[:,0], yvals[:,0], zvals[:,0], marker='.', c='k', alpha=0.15, s=10) # Initial posits
        ax.scatter(xend, yend, zend, marker='.', c='k', alpha=0.15, s=10) # Final posits
        ax.scatter(xvals[i,0], yvals[i,0], zvals[i,0], marker='.', c='k', alpha=0.5, s=10) # Select Initial posits
        ax.scatter(xend[i], yend[i], zend[i], marker='.', c='k', alpha=0.5, s=10) # Select Final posits
        ax.set_aspect('equal')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.elev=4
        ax.azim=25
        ax.set_title("3D view of 1 particle")
    
#        ax2.plot(posxyz[1,i,:].T, np.rad2deg(theta[:10,:].T))
#        ax2.set_title("Theta vs. Y depth (1 particle)")
#        ax2.set_ylabel("Theta (degs)")
#        ax2.set_xlabel("Spatial Y")

        ax3.clear()
        ax3.plot(dtheta[i,:].T, c='r')
        ax3.set_ylim(-0.006, 0.006)
        ax3.set_title("dTheta vs. Step # (1 particle) (TOFAdj = " + str(tofadj[i]) + ")")
        ax3.set_ylabel("dTheta (drad)")
    
#        ax4.plot(pathl[i,:].T, np.rad2deg(theta[:10,:].T))
#        ax4.set_title("Theta vs. Path length (1 particle)")
#        ax4.set_ylabel("Theta (degs)")
#        ax4.set_xlabel("Path length")
        fig.savefig(os.path.join(r'C:\Users\Scott\Documents\TDYNO\plots\2016-12-16 pmov', str(i).zfill(5) + '.png'))

