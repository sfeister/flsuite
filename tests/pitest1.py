#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
pitest1.py: Test of reading a proton imaging plotfile contents

Created by Scott Feister on Wed Dec 07 15:31:48 2016

Changelog:
    2016-12-12 Spun off this program from trimtest.py (analysis of TRIM data)
"""

import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":
    fn = r"C:\Users\Scott\Documents\temp\subdata.hdf5"
    f = h5py.File(fn, 'r')
    ix_arr = np.argsort(f["ProtonData"][:,0], kind='mergesort') # Mergesort keeps the paths in order, for a given ID # # TODO: Sort in batch sizes?
    dat = f["ProtonData"][...]
    dat = dat[ix_arr, :] # Sorted and transposed

    
    # http://stackoverflow.com/questions/19125661/find-index-where-elements-change-value-numpy
    idarr = dat[:, 0] # Array of particle ID numbers
    splitix = np.where(idarr[:-1] != idarr[1:])[0] + 1 # Array of indices at which particle ID changes (split here)
    edgeix = np.hstack((0, splitix, len(idarr))) # Array of start and stop indices; N + 1 long, where N is number of particle tracks
    traklens = np.diff(edgeix) # Array of particle track lengths
    ntraks = len(traklens) # Number of particle tracks
    tracks = np.split(dat, splitix, axis=0) # Original array, split up into discrete lists by particle
    #plt.figure(1)
    #plt.clf()
    tlmax = 500 # Max track length; Either the max number in the arrays, or 50 steps
    trakarr = np.empty((ntraks, tlmax, dat.shape[1]))
    trakarr.fill(np.nan)
    
    for i in range(ntraks):
        trakarr[i,:min(traklens[i], tlmax),:] = dat[edgeix[i]:edgeix[i] + min(traklens[i], tlmax),:]
    
    dtrakarr = np.diff(trakarr, axis=1) # Difference along the traveling axis
    dr = np.sqrt(dtrakarr[:,:,1]**2 + dtrakarr[:,:,2]**2 + dtrakarr[:,:,3]**2) # N x Steps-1 long, Path length change between steps
    pathl = np.concatenate((np.zeros((dr.shape[0],1)), np.cumsum(dr, axis=1)), axis=1) # N x Steps long, Total path length up to this step

    fig = plt.figure(1)
    fig.clear()
    ax = fig.add_subplot(111)

    ntplot = 400 # Number of tracks to plot
    maxpathl = 0.6
    
    ct = pathl[:ntplot,:] > maxpathl
    xvals = trakarr[:ntplot,:,1]
    yvals = trakarr[:ntplot,:,2]
    zvals = trakarr[:ntplot,:,3]
    xvals[ct] = np.nan
    yvals[ct] = np.nan
    zvals[ct] = np.nan
    
    fig = plt.figure(1)
    fig.clear()
    ax = fig.add_subplot(111)
    ax.plot(xvals.T, zvals.T)
    ax.set_aspect('equal')
    #ax.set_ylim(-3000, 3000)
    #ax.set_xlim(0.27, 0.29)
    """
    fig = plt.figure(1)
    fig.clear()
    ax = fig.add_subplot(111)
    fig2 = plt.figure(3)
    fig2.clear()
    ax2 = fig2.add_subplot(111, projection='3d')

    for track in tracks:
        ax.plot(track[:,2], track[:,3])
        ax2.plot(track[:,2], track[:,3], track[:,4], c='k', alpha=0.01)
    #trackarr = np.stack(tracks)

    ax2.set_ylim(-3000, 3000)
    ax2.set_zlim(-3000, 3000)

    fig = plt.figure(2)
    fig.clear()
    ax2 = fig.add_subplot(111, projection='3d')
    ax2.scatter(dat[:,2], dat[:,3], dat[:,4], marker='.', alpha=0.05, s=1)
    ax2.set_ylim(-3000, 3000)
    ax2.set_zlim(-3000, 3000)
    """