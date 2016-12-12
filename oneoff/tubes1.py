#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
tubes1.py: Attempt to reconstruct particle diffusion info from vtk array

Created by Scott Feister on Wed Dec 07 11:46:50 2016
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":
    zvals = np.arange(0, 100)
    xvals = zvals**2
    yvals = zvals**1.5
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xvals, yvals, zvals)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    
