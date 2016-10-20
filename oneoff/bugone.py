#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
bugone.py: Load a FLASH file and debug it

Created by Scott Feister on Wed Oct 19 18:58:45 2016

USAGE:
python -i bugone.py     (Loads and gives you ds to play with)
"""

import yt
from flsuite.flyt import tstools as tst
import flsuite.sftools as sf
from flsuite.flyt import flyt
#import flsuite.flyt.morefields # Add a few extra, custom derived fields to FLASH datasets loaded through YT
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime

if __name__ == "__main__":
    # Note that this will run in parallel, so no prints please
    datdir = r'/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND' # Directory holding HDF5 FLASH outputs
    basenm = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
    fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_????') # yt-time-series filename pattern for plot files

    outroot = r'/home/sfeister/myouts' # Folder to store outputs
    outdir = sf.subdir(outroot,"scratch_dev")

    # Three objects you may want to use
    ts = yt.load(fnpatt)
    ds = ts[100]
    dd = ds.all_data()
    
    