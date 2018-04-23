#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parLaser.py: Defines a class for writing laser strings into the FLASH file

Created by Scott Feister on Wed Feb 14 13:39:38 2018

TODO: Add further documentation. -SKF 4-23-2018
"""

import numpy as np

class parLaser:
    """ A class containing the parameters of a single beam and pulse flash.par (runtime parameters) input """
    
    def __init__(self, lasnum, laslbl=None):
        # Special variables for this class
        self.lasnum = int(lasnum) # A number for the laser, considered both the beam-number and pulse-number
        self.laslbl = laslbl # A label for the laser, e.g. 'Quad24', which is put into the title
        self.lens = None # A 3-element list or array with values for [lensx, lensy, lensz]
        self.targ = None # A 3-element list or array with values for [targetx, targety, targetz]
        self.times = None # Will become an array of times for the laser pulse
        self.powers = None # Will become an array of powers for the laser pulse

        # Basically initializes everything else found under the "ed_XXXX_1" header
        self.gridNAngularTics = None
        self.gridNSemiAxisMajorTics = None
        self.gridNSemiAxisMinorTics = None
        self.numberOfRays = None
        self.gaussianCenterMajor = None
        self.gaussianCenterMinor = None
        self.gaussianExponent = None
        self.gaussianRadiusMajor = None
        self.gaussianRadiusMinor = None
        self.gridDeltaSemiAxisMajor = None
        self.gridDeltaSemiAxisMinor = None
        self.initialRaySpeed = None
        self.lensSemiAxisMajor = None
        self.wavelength = None
        self.semiAxisMajorTorsionAngle = None
        self.targetSemiAxisMajor = None
        self.targetSemiAxisMinor = None
        self.crossSectionFunctionType = None
        self.gridType = None
        self.semiAxisMajorTorsionAxis = None
        self.ignoreBoundaryCondition = None
    
    def makePar(self):
        par = ''
        # Check for valid laser number (>0, integer for FLASH)
        if self.lasnum < 1:
            raise Exception('Cannot write pulse or beam: Laser number "lasnum" must be an integer greater than zero.')
            
        ## WRITE THE PULSE
        # Some quick checks
        if not (hasattr(self.powers, '__len__') and (not isinstance(self.powers, str))):
            raise Exception('Cannot write pulse: Powers are not specified as a list or 1D array.')
        if not (hasattr(self.times, '__len__') and (not isinstance(self.times, str))):
            raise Exception('Cannot write pulse: Times are not specified as a list or 1D array.')
        if len(self.powers) != len(self.times):
            raise Exception("Cannot write pulse {}: Powers and times have different numbers of elements.".format(self.lasnum))
        
        # Write the pulse header and values
        par += "\n"
        if self.laslbl is not None:
            par += "## Begin pulse {} ({}):\n".format(self.lasnum, self.laslbl)
        else:
            par += "## Begin pulse {}:\n".format(self.lasnum)

        par += "ed_numberOfSections_{} = {}\n".format(self.lasnum, len(self.powers))
        for i, power in enumerate(self.powers, start=1):
            par += "ed_power_{}_{} = {}\n".format(self.lasnum, i, power)
        for i, time in enumerate(self.times, start=1):
            par += "ed_time_{}_{} = {}\n".format(self.lasnum, i, time)
                
        ## WRITE THE BEAM
        # Do some checks for lens and target arrays
        if not (hasattr(self.lens, '__len__') and (not isinstance(self.lens, str))):
            raise Exception('Cannot write beam: Lens is not specified as a 3-element list or array.')
        if not (hasattr(self.targ, '__len__') and (not isinstance(self.targ, str))):
            raise Exception('Cannot write beam: Targ is not specified as a 3-element list or array.')
        if len(self.lens) != 3:
            raise Exception('Cannot write beam: Lens has less or more than 3 elements.')
        if len(self.targ) != 3:
            raise Exception('Cannot write beam: Targ has less or more than 3 elements.')
        
        # Write the beam (lens and target)
        par += "\n"
        if self.laslbl is not None:
            par += "## Begin beam {} ({}):\n".format(self.lasnum, self.laslbl)
        else:
            par += "## Begin beam {}:\n".format(self.lasnum)
        
        par += "ed_pulseNumber_{} = {}\n".format(self.lasnum, self.lasnum)
        
        for i, dim in enumerate(["X", "Y", "Z"]):
            par += "ed_lens{}_{} = {}\n".format(dim, self.lasnum, self.lens[i])
        for i, dim in enumerate(["X", "Y", "Z"]):
            par += "ed_target{}_{} = {}\n".format(dim, self.lasnum, self.targ[i])
    
        # Write any remaining beam items (anything not set to 'None')
        goodkeys = set(self.__dict__.keys()) - set(["lasnum", "laslbl", "lens", "targ", "powers", "times"])
        for key in goodkeys:
            if getattr(self, key) is not None:
                par += "ed_{}_{} = {}\n".format(key, self.lasnum, getattr(self, key))
        
        return par
    
    def __str__(self):
        return self.makePar()
        
if __name__ == "__main__":
    ## EXAMPLE USAGE
    
    las1 = parLaser(1, laslbl="Red laser")
    las1.lens = [20, 20, 30]
    las1.targ = [20, 30, 40]
    las1.powers = np.array([1,2,3,4,5])
    las1.times = np.array([10,11,12,13,14])
    las1.wavelength = 0.808
    
    las2 = parLaser(2, laslbl="Many-rays Laser")
    las2.lens = [15, 15, 23]
    las2.targ = [22, 22, 41]
    las2.powers = np.array([1,2,3,4,5,6,7])
    las2.times = np.array([10,11,12,13,14,15,16])
    las2.numberOfRays = 10000
    
    las3 = parLaser(3, laslbl="Gaussian Laser")
    las3.lens = [14, 14, 16]
    las3.targ = [40, 50, 52]
    las3.powers = np.array([2,2.5,3])
    las3.times = np.array([10,11,12])
    
    las3.crossSectionFunctionType = "gaussian2D" # 2D Gaussian Beam
    las3.gaussianExponent = 4.0 # 4.0 for supergaussian profile
    las3.gaussianRadiusMajor = 0.048
    las3.gaussianRadiusMinor = 0.048

    par = las1.makePar()
    par += las2.makePar()
    par += las3.makePar()

    print(par)
    