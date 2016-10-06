#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
tstools.py: Tools for manipulating YT Time Series (bundles of FLASH datasets)

Created by Scott Feister on Wed Oct 05 18:18:14 2016

Example usage:
import flsuite.flyt.tstools as tst

CHANGELOG:
2016-10-05 Created tstools.py, filled with some functions from tstest3.py
"""

import numpy as np

# TODO: Parallelize this code
def tsinf(ts):
	""" Walks through the time series and returns a dict with some basic info like time steps """
	numfiles = len(ts)

	ts_inf = {}
	ts_inf['times_ns'] = np.zeros(numfiles)
	
	for i in range(numfiles):
		ts_inf['times_ns'][i] = ts[i].current_time.in_units('ns').v

	return ts_inf


if __name__ == "__main__":
    pass
