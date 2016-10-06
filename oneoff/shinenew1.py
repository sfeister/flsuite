#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
shinenew1.py: Do shinethru analysis for TDYNO sim, but use flsuite tools

Created by Scott Feister on Wed Oct 05 18:01:01 2016

CHANGELOG:
2016-10-05 Created shinenew1.py. Copied many things from shinethru3.py and tstest3.py.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg') # Headless plotting
import matplotlib.pyplot as plt
import yt
yt.enable_parallelism() # Tap into yt's mpi4py parallelism (e.g. now can call via mpirun -np 10 python <blah>.py)
yt.funcs.mylog.setLevel(30) # This sets the output notification threshold to 30, WARNING. Default is 20, INFO.
import flsuite as fl
import flsuite.flyt.morefields # Add a few extra, custom derived fields to FLASH datasets loaded through YT

#TODO: Make this modular, so you can get what you need out of it. Could it return a callback?
def anlzD(ds):
	"""Perform custom analysis of a single yt dataset
	
	Returns the analysis data packaged up into a dictionary, ds_anls
	"""
	ds_anls = {}
	ds_anls['times_ns'] = ds.current_time.in_units('ns').v

	reg = ds.all_data() # Region of interest
	ds_anls['Bav_code'] = reg.mean("magtot").v
	
	return ds_anls

#TODO: This will eventually be made modular, but tied in with anlzD
def plotT(ts_anls, plotdir='.'):
	"""Generate and save custom plot(s) of the results of a anlzT() time series analysis
	
	ts_anls: Output of anlzT()
	"""
	plt.figure(1)
	plt.clf()
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_code']*np.sqrt(4*np.pi)*1e-3, label='All space')
	plt.xlabel("Simulation time (ns)")
	plt.ylabel("Mean magnetic field strength (kG)")
	plt.title("TDYNO_BAND, Magnetic field through time")
	plt.legend(loc='upper left')
	plt.savefig(os.path.join(plotdir,'shinethru1.png'), dpi=300)

	return 0 # TODO: Return the figure handles?

# TODO: Make this generic so that it works regardless of anlzD and plotT
def anlzT(ts):
	""" Walks through the time series and does a custom analysis "anlzone(ds)" IN PARALLEL on all datasets in time series
	Follows walkthrough at: http://yt-project.org/doc/analyzing/parallel_computation.html?highlight=parallel#parallelization-over-multiple-datasets-including-time-series
	"""
	numfiles = len(ts)

	if yt.is_root():
		print "Number of files to analyze:", numfiles
		
	# Define an empty storage dictionary for collecting information
	# in parallel through processing
	storage = {}

	for sto, ds in ts.piter(storage=storage):
		print "Reading file:", ds
		ts_anls_mini = anlzD(ds)
		
		sto.result = ts_anls_mini
		sto.result_id = str(ds)
	
	# Repackage in serial (could be better done)
	ts_anls = {} # Make sure everything above gets defined here!
	ts_anls['times_ns'] = np.zeros(numfiles)
	ts_anls['Bav_code'] = np.zeros(numfiles)

	print storage
	
	sortd = sorted(storage.keys()) # Sorted result_id keys
	for i in range(numfiles):
		ts_anls_mini = storage[sortd[i]]
		for k in ts_anls_mini:
			ts_anls[k][i] = ts_anls_mini[k]

	return ts_anls


if __name__ == "__main__":
	# Note that this will run in parallel, so no prints please
	datdir = r'/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND' # Directory holding HDF5 FLASH outputs
	basenm = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
	fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_??[0,5]0') # yt-time-series filename pattern for plot files

	plotdir = r'/home/sfeister/myouts/TDYNO_BAND'
	ts = yt.load(fnpatt) # Load a time series
	ts_anls = anlzT(ts)
	
	if yt.is_root():
		print "Plotting..."
		print ts_anls
		plotT(ts_anls, plotdir=plotdir)
