#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
laserexamples.py: Show two different ways to create laser strings for a flash.par file

Created by Scott Feister on Fri Jul 27 15:37:49 2018
"""
  
from flsuite.parLaser import parLaser, parLasers
import numpy as np

# Example 1: Three lasers, each with different parameters
las1 = parLaser(1, laslbl="Second Laser, 808 nm")
las1.lens = [20, 20, 30]
las1.targ = [20, 30, 40]
las1.powers = np.array([1,2,3,4,5])
las1.times = np.array([10,11,12,13,14])
las1.wavelength = 0.808

las2 = parLaser(2, laslbl="Second Laser, Many rays")
las2.lens = [15, 15, 23]
las2.targ = [22, 22, 41]
las2.powers = np.array([1,2,3,4,5,6,7])
las2.times = np.array([10,11,12,13,14,15,16])
las2.numberOfRays = 10000

las3 = parLaser(3, laslbl="Third Laser, Gaussian profile")
las3.lens = [14, 14, 16]
las3.targ = [40, 50, 52]
las3.powers = np.array([2,2.5,3])
las3.times = np.array([10,11,12])

las3.crossSectionFunctionType = "gaussian2D" # 2D Gaussian Beam
las3.gaussianExponent = 4.0 # 4.0 for supergaussian profile
las3.gaussianRadiusMajor = 0.048
las3.gaussianRadiusMinor = 0.048

las1.write('laser1.txt', 'w')
las2.write('laser2.txt', 'w')
las3.write('laser3.txt', 'w')
print("Three lasers written to 'laser1.txt', 'laser2.txt', 'laser3.txt'")

## Example 2: Make a whole bunch of the same laser, but with lens at varying x value
# Uses the "parlasers" class
l = parLasers(10) # Laser list

for i in range(len(l)):
    l[i].lens = [i*10, 0, 0] # This is the only thing changing between the ten lasers!
    l[i].targ = [5, 5, 5]
    l[i].powers = np.array([1,2,3,4,5])
    l[i].times = np.array([10,11,12,13,14])
    l[i].numberOfRays = 10000
    l[i].crossSectionFunctionType = "gaussian2D" # 2D Gaussian Beam
    l[i].gaussianExponent = 4.0 # 4.0 for supergaussian profile
    l[i].gaussianRadiusMajor = 0.048
    l[i].gaussianRadiusMinor = 0.048

l.write('tenlasers.txt', 'w')
print("Ten lasers written to 'tenlasers.txt'")
