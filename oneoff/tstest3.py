#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
tstest3.py: Time-series test for yt. Analyzing B fields at core of TDYNO

Created by Scott Feister around Oct 1 2016

Example usage:
mpirun -np 16 python tstest2.py

CHANGELOG:
2016-10-06 Copied contents of tstest2.py on Mira into this folder
"""

# Version 2 tests out parallel processing using mpi4py
# 

print "Importing yt..."
import yt
yt.enable_parallelism()
#from yt import YTArray
print "Performing other imports..."
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import scottderived # Local script adds magnetic and velocity magnitudes to loaded fields

#print "Beginning script."

def tsinfo(ts):
	""" Walks through the time series and returns a dict with some basic info like time steps """
	numfiles = len(ts)

	ts_inf = {}
	ts_inf['times_ns'] = np.zeros(numfiles)
	
	for i in range(numfiles):
		ts_inf['times_ns'][i] = ts[i].current_time.in_units('ns').v

	return ts_inf

def tsanalyze(ts):
	""" Walks through the time series and does a custom analysis """
	numfiles = len(ts)

	ts_anls = {}
	ts_anls['times_ns'] = np.zeros(numfiles)
	ts_anls['Bav_code'] = np.zeros(numfiles)
	
	for i in range(numfiles):
		ds = ts[i]
		print "Reading file:", ds
		reg = ds.all_data() # Region of interest #TODO: Define a disk, sphere, etc.
		ts_anls['times_ns'][i] = ds.current_time.in_units('ns').v
		ts_anls['Bav_code'][i] = reg.mean("magtot").v

	return ts_anls

def tsanpar(ts):
	""" Walks through the time series and does a custom analysis IN PARALLEL
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
		ts_anls_mini = {}
		ts_anls_mini['times_ns'] = ds.current_time.in_units('ns').v

		reg = ds.all_data() # Region of interest #TODO: Define a disk, sphere, etc.
		ts_anls_mini['Bav_code'] = reg.mean("magtot").v

		reg2 = ds.sphere(ds.domain_center, (500, 'um')) # 500-um-radius sphere about the center
		ts_anls_mini['Bav_cent_code'] = reg2.mean("magtot").v

		reg3 = ds.sphere(ds.domain_center, (100, 'um')) # 100-um-radius sphere about the center
		ts_anls_mini['Bav_cent2_code'] = reg3.mean("magtot").v

		reg5 = ds.disk(ds.arr([0,-0.3125,0.4],'cm'),[0,1,0],ds.quan(.018,'cm'),ds.quan(0.3125*2,'cm')) # 180 um radius cylinder, full-length of domain (1 deg aperture equivalent)
		ts_anls_mini['Bav_cyl1_code'] = reg5.mean("magtot").v

		reg6 = ds.disk(ds.arr([0,-0.3125,0.4],'cm'),[0,1,0],ds.quan(.055,'cm'),ds.quan(0.3125*2,'cm')) # 550 um radius cylinder, full-length of domain (3 deg aperture equivalent)
		ts_anls_mini['Bav_cyl2_code'] = reg6.mean("magtot").v

		reg7 = ds.disk(ds.arr([0,-0.3125,0.4],'cm'),[0,1,0],ds.quan(.091,'cm'),ds.quan(0.3125*2,'cm')) # 910 um radius cylinder, full-length of domain (5 deg aperture equivalent)
		ts_anls_mini['Bav_cyl3_code'] = reg7.mean("magtot").v

		reg8 = ds.disk(ds.arr([0,-0.3125,0.4],'cm'),[0,1,0],ds.quan(.182,'cm'),ds.quan(0.3125*2,'cm')) # 1.82 mm radius cylinder, full-length of domain (CR-39 equivalent)
		ts_anls_mini['Bav_cyl4_code'] = reg8.mean("magtot").v
		
		reg9 = ds.disk(ds.arr([0,-0.3125,0.4],'cm'),[0,1,0],ds.quan(.346,'cm'),ds.quan(0.3125*2,'cm')) # 3.46 mm radius cylinder, full-length of domain (CR-39 equivalent)
		ts_anls_mini['Bav_cyl5_code'] = reg9.mean("magtot").v
	

	
		proj = yt.ProjectionPlot(ds, "x", ["bdry", "magtot"], data_source=reg7)
		proj.save("vids")

		proj = yt.ProjectionPlot(ds, "y", ["bdry", "magtot"], data_source=reg7)
		proj.save("vids")

		sto.result = ts_anls_mini
		sto.result_id = str(ds)
	
	# Repackage in serial (could be better done)
	ts_anls = {} # Make sure everything above gets defined here!
	ts_anls['times_ns'] = np.zeros(numfiles)
	ts_anls['Bav_code'] = np.zeros(numfiles)
	ts_anls['Bav_cent_code'] = np.zeros(numfiles)
	ts_anls['Bav_cent2_code'] = np.zeros(numfiles)
	ts_anls['Bav_cyl1_code'] = np.zeros(numfiles)
	ts_anls['Bav_cyl2_code'] = np.zeros(numfiles)
	ts_anls['Bav_cyl3_code'] = np.zeros(numfiles)
	ts_anls['Bav_cyl4_code'] = np.zeros(numfiles)
	ts_anls['Bav_cyl5_code'] = np.zeros(numfiles)

	print storage
	
	sortd = sorted(storage.keys()) # Sorted result_id keys
	for i in range(numfiles):
		ts_anls_mini = storage[sortd[i]]
		for k in ts_anls_mini:
			ts_anls[k][i] = ts_anls_mini[k]

	return ts_anls

# Suppress 'yt: [Info]' outputs to terminal"
yt.funcs.mylog.setLevel(30) # This sets the output notification threshold to 30, WARNING. Default is 20, INFO.

# Load FLASH plotfiles as yt time-series
datdir = r'/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND' # Directory holding HDF5 FLASH outputs
basenm = r'tdyno2016_' # Prefix for plot filenames, which is defined as "basenm" in the flash.par file
fnpatt = os.path.join(datdir, basenm + 'hdf5_plt_cnt_???[0,2,4,6,8]') # yt-time-series filename pattern for plot files
# hdf5_plt_cnt_???0
# hdf5_plt_cnt_??00
# hdf5_plt_cnt_???[0,2,4,6,8]
# hdf5_plt_cnt_????

ts = yt.load(fnpatt) # Load FLASH outputs as yt time-series (Doesn't actually load data into RAM; just creates references)

#ts_inf = tsinfo(ts)

ts_anls = tsanpar(ts)
if yt.is_root():
	print "I am root. Here is the output."
	print ts_anls
	
	print "Plotting..."
	plt.figure(1)
	plt.clf()
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_code']*np.sqrt(4*np.pi)*1e-3, label='All space')
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_cent_code']*np.sqrt(4*np.pi)*1e-3, label='Central sphere, r=0.5 mm')
	#plt.plot(ts_anls['times_ns'], ts_anls['Bav_cent2_code']*np.sqrt(4*np.pi)*1e-3, label='Central sphere, r=0.1 mm')
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_cyl1_code']*np.sqrt(4*np.pi)*1e-3, label=r'Cyl., r=0.18 mm ($2\theta=1^\circ$)')
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_cyl2_code']*np.sqrt(4*np.pi)*1e-3, label=r'Cyl., r=0.55 mm ($2\theta=3^\circ$)')
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_cyl3_code']*np.sqrt(4*np.pi)*1e-3, label=r'Cyl., r=0.91 mm ($2\theta=5^\circ$)')
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_cyl4_code']*np.sqrt(4*np.pi)*1e-3, label=r'Cyl., r=1.82 mm ($2\theta=10^\circ$)')
	plt.plot(ts_anls['times_ns'], ts_anls['Bav_cyl5_code']*np.sqrt(4*np.pi)*1e-3, label=r'Cyl., r=3.46 mm ($2\theta=19^\circ$)')
	plt.xlabel("Simulation time (ns)")
	plt.ylabel("Mean magnetic field strength (kG)")
	plt.title("TDYNO_BAND, Magnetic field through time")
	plt.legend(loc='upper left')
	plt.savefig('btime3_2.png', dpi=300)
	
#ts_anls = tsanalyze(ts)