#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
oct16pi1.py: Custom analysis for proton imaging, oct 2016

Created by Scott Feister on Mon Oct 17 12:02:21 2016

Analyze the SCRIPT1/RUN1 and SCRIPT1/RUN2 directories
"""

import os
from flsuite.PI import piHugeAnalysis
import sftools as sf

if __name__ == "__main__":
    PIroot = r"/home/sfeister/tzef/NIF/TDYNO_BAND/PINHOLE/SCRIPT1"
    outroot = r"/gpfs/mira-home/sfeister/myouts/TDYNO_BAND/PI"
    dirnames = [r"RUN1", r"RUN2"]
    simnames = [r"TDYNO_BAND_elec", r"TDYNO_BAND_no_elec"]
    for i in range(len(dirnames)):
        PIdir = os.path.join(PIroot, dirnames[i])
        outdir = sf.subdir(outroot, simnames[i])
        print("Analyzing directory: " + PIdir)
        piHugeAnalysis(PIdir, basenm=r"tdyno2016_", simname=simnames[i], outdir=outdir)
