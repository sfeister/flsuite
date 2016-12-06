#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
PIhdftest.py: Test read of new proton imaging HDF5 output (in the output plotfile)

Created by Scott Feister on Tue Dec 06 16:29:48 2016
"""
import os
import matplotlib as mpl
mpl.use("Agg")
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":
    fn = r"/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_PROTONS/133kJ/20ns/3MeV/tdyno2016_hdf5_plt_cnt_0001"
    #fn = r"/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_PROTONS/133kJ/20ns/tdyno2016_hdf5_chk_0320"
    f = h5py.File(fn, 'r')
    for k in f.keys():
        print (str(k) + ": " + str(f[k]))
        
    print("Extracting data...")
    #f["ProtonData"][:,:]
    outdir = r"/home/sfeister/myouts/TDYNO_PROTONS"
    fig = plt.figure(1)
    plt.clf()
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.hist(f['ProtonData'][:,i], 50)
        plt.title("Index " + str(i))
    plt.suptitle("ProtonData Histograms")
    fig.savefig(os.path.join(outdir, "myplot.png"), dpi=150)

    print("Making second plot...")
    fig = plt.figure(2)
    plt.clf()
    for i in range(4):
        plt.subplot(2,2,i+1)
        plt.plot(f['ProtonData'][:,i])
        plt.title("Index " + str(i))
        #plt.xlim(50000, 52000)
    plt.suptitle("ProtonData plotted in sequence")
    fig.savefig(os.path.join(outdir, "myplot2.png"), dpi=150)
    
    print("Making 3D plot...")
    # 3D Scatter plot
    fig = plt.figure(3)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(f['ProtonData'][::500,1],f['ProtonData'][::500,2],f['ProtonData'][::500,3], s=10)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.savefig(os.path.join(outdir, "myplot3.png"), dpi=300)
