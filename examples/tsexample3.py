#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tsexample3.py: Example, Make customized matplotlib density slice plot at every timestep (a png movie)

Modify the file pattern for your HDF5 files ("filepatt" variable below), then run via "python tsexample3.py".
If mpi4py is installed on your system, you could run this in parallel on four cores by "mpirun -np 4 python tsexample3.py"

Created by Scott Feister on Fri Jul 27 09:38:35 2018
"""

import os, sys
import numpy as np
import yt
import flsuite.flyt.tstools as tst
from flsuite.flyt.flyt import get_simdata
import matplotlib.pyplot as plt
from matplotlib import colors

def anlzD(ds, anlsD, outdir):
    """ Use yt to make plots (z-axis slice) of density and electron temperature for this timestep """
    pltnum = str(ds)[-4:] # Plot number string, e.g. strips "0020" off of "lasslab_hdf5_plt_cnt_0020"

    # Extract underlying grid data into uniformly spaced arrays
    cg = get_simdata(ds)
    [xmin, ymin] = cg.left_edge[:2].v
    [xmax, ymax] = cg.right_edge[:2].v
    Dens = np.squeeze(cg['dens'].v) # 2D array of cell density values

    # Make a custom plot of density using matplotlib
    fig, ax = plt.subplots()
    vmin = 1e-4
    vmax = 1e-1
    norm = colors.LogNorm(vmin, vmax) # Log scale
    #norm = colors.Normalize(vmin, vmax) # Linear scale
    im = ax.pcolorfast([xmin, xmax], [ymin, ymax], Dens.T, norm=norm, cmap='plasma')
    ax.set_title("My custom plot of Density")
    ax.set_aspect("equal")
    ax.set_xlabel("X (cm)")
    ax.set_ylabel("Y (cm)")
    plt.colorbar(im, label="Density (g/cc)")
    
    fig.savefig("dens_" + pltnum + ".png")
    plt.close(fig)
    
if __name__ == "__main__":
    filepatt = '/home/scott/myouts/lasslab_hdf5_plt_cnt_????' # File pattern (note the ????)
    ts = yt.load(filepatt)
    anlsT = tst.anlzT(ts, anlzD, outdir='.')
