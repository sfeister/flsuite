#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
nifpi.py: Energy shift comparison for proton imaging

Created by Scott Feister on Wed Dec 14 10:50:50 2016
"""

import flsuite.PI as PI

if __name__ == "__main__":
    PIdir1 = r"/projects/CosmicLaser/tzeferac/NIF/TDYNO_PROTONS/PR_SPECTRUM/133 kJ"
    PIdir2 = r"/projects/CosmicLaser/tzeferac/NIF/TDYNO_PROTONS/PR_SPECTRUM/200 kJ"
    outdir1 = r"/home/sfeister/myouts/NIF_comparisons/ProtonSpectrum/133 kJ 22_5 ns"
    outdir2 = r"/home/sfeister/myouts/NIF_comparisons/ProtonSpectrum/200 kJ 20_0 ns"

    PI.piHugeAnalysis(PIdir1, basenm=r"tdyno2016_", simname='NIF 133 kJ', outdir=outdir1)
    PI.piHugeAnalysis(PIdir2, basenm=r"tdyno2016_", simname='NIF 200 kJ', outdir=outdir2)