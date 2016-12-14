#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
PITest1.py: Get the proton radiography working for "electron spectrum" style outputs

Created by Scott Feister on Fri Oct 21 12:12:56 2016
"""

import flsuite
import sftools as sf
import numpy as np

if __name__ == "__main__":
    PIdir = r"C:\Users\Scott\Documents\temp\oct2016\PI 200kjSelf"
    basenm = r"tdyno2016PI_"
    fns = sf.getfns(PIdir, prefix = basenm + 'ElectronSpectrometryDetector')
    Es = np.genfromtxt
    pass
