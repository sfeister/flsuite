#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
tstools.py: Tools for manipulating YT Time Series (bundles of FLASH datasets)

Created by Scott Feister on Wed Oct 05 18:18:14 2016

Example usage:
import flsuite.flyt.tstools as tst

CHANGELOG:
2016-10-05 Created tstools.py, filled with some functions from tstest3.py
2016-10-13 Greatly expanded the scope, automating post-analysis.

TODO:
* Tie together plotT and anlsD, perhaps making a class!
* Enable further modularity to allow multiple instances of anlsD, plotT to be included (e.g. allows a list of functions for each)
   Possibly, even have the anlsD allow anlsT inputs (for auto-scaling, etc.). Or, allow classes manually setting parameters like mins, maxes, etc. Could become a monstrosity if that, though!
"""

import os
import numpy as np
import cPickle as pickle
import matplotlib
matplotlib.use('Agg') # Headless plotting
import matplotlib.pyplot as plt
import yt
yt.enable_parallelism() # Tap into yt's mpi4py parallelism (e.g. now can call via mpirun -np 10 python <blah>.py)
yt.funcs.mylog.setLevel(30) # This sets the output notification threshold to 30, WARNING. Default is 20, INFO.
import flsuite as fl
import flsuite.flyt.morefields # Add a few extra, custom derived fields to FLASH datasets loaded through YT
import sftools as sf
import numbers

# TODO: Incorporate this into the new analysis schema
def tsinf(ts):
	""" Walks through the time series and returns a dict with some basic info like time steps """
	numfiles = len(ts)

	ts_inf = {}
	ts_inf['times_ns'] = np.zeros(numfiles)
	
	for i in range(numfiles):
		ts_inf['times_ns'][i] = ts[i].current_time.in_units('ns').v

	return ts_inf


def example_anlzD(ds, outdir='.'):
    """(*AnaLyZe Dataset*) Perform custom analysis of a single yt dataset 
    
    Returns the analysis data packaged up into a dictionary, anlsD
    """
    anlsD = {}
    anlsD['times_ns'] = ds.current_time.in_units('ns').v

    reg = ds.all_data() # Region of interest
    anlsD['Bav_code'] = reg.mean("magtot").v
    
    # Slice plots
    slc = yt.SlicePlot(ds, 'x', ['density', 'velz', 'magtot'])
    slc.set_log('velz', False)
    slc.set_cmap('velz', 'RdBu')
    slc.set_unit('velz', 'um/ps')
    slc.set_zlim('velz', -1.0, 1.0)
    slc.set_zlim('density', 1e-6, 1e1)
    slc.set_font({'family': 'sans-serif', 'style': 'italic',
        'weight': 'bold', 'size': 24})
    # Projection plots
    slc.annotate_title('Slice at X=0')
    slc.save(sf.subdir(outdir, 'xslice')) # Known issue: Technically, this can be a race condition among processors of creating this sub-directory. In practice, it's probably fine most of the time.

    proj = yt.ProjectionPlot(ds, 'x', ['depo'], method='mip')
    proj.set_zlim('depo', 1e-15, 1e14)
    proj.set_font({'family': 'sans-serif', 'style': 'italic',
        'weight': 'bold', 'size': 24})
    proj.annotate_title('Projection through X, maximum value')
    proj.save(sf.subdir(outdir, 'xmip'))

    return anlsD

def example_plotT(anlsT, plotdir='.'):
    """(*Plot Time series*) Generate and save custom plot(s) of the results of a anlzT() time series analysis
    
    anlsT: Output of anlzT()
    """
    plt.figure(1)
    plt.clf()
    plt.plot(anlsT['times_ns'], anlsT['Bav_code']*np.sqrt(4*np.pi)*1e-3, label='All space')
    plt.xlabel("Simulation time (ns)")
    plt.ylabel("Mean magnetic field strength (kG)")
    plt.title("TDYNO_BAND, Magnetic field through time")
    plt.legend(loc='upper left')
    plt.savefig(os.path.join(plotdir,'shinethru1.png'), dpi=300)

    return 0 # TODO: Return the figure handles?

# TODO: Homogenize string outputs from anlsD such that they become numpy arrays (?)
def anlzT(ts, anlzD, outdir='.', plotT=None):
    """(*AnaLyZe Time series*) Walks through the time series and does a custom analysis "anlzone(ds)" IN PARALLEL on all datasets in time series
    Follows walkthrough at: http://yt-project.org/doc/analyzing/parallel_computation.html?highlight=parallel#parallelization-over-multiple-datasets-including-time-series
    
    anlzD: Function defined with inputs (ds, outdir='.') and output anlsD (dict of values with keys). Run on every dataset.
    plotT: (Optional) Function defined with inputs (anlsT, outdir='.') and no particular output. Run once at end of analysis.
    
    anlsD: Dictionary compiling output from a single dataset. Output of anlzD. 
    anlsT: Dictionary compiling outputs from time-series, keys matching those of anlsD.
    """
    numfiles = len(ts)

    if yt.is_root():
        print "Number of files to analyze:", numfiles
        
    ## Analyze datasets using yt parallelization
    storage = {} # Define an empty storage dictionary for collecting information in parallel through processing

    for sto, ds in ts.piter(storage=storage): # Iterate over datasets ds in time-series ts
        print "Reading file:", ds
        anlsD = anlzD(ds, outdir=outdir) # It would be nice not to require the outdir argument
        
        if not isinstance(anlsD, dict): # If for some reason it returned something other than a dictionary
            anlsD = {'anlzD_out':anlsD}
        
        sto.result = anlsD
        sto.result_id = str(ds)
        
    sortd = sorted(storage.keys()) # Sorted result_id keys

    ## Re-organize results in serial
    # Pre-allocate output dictionary anlsT, using keys from anlsD
    # Intelligently pre-allocate a time-series analysis dictionary (anlsT) using keys and values of the first single-dataset analysis dictionary (anlsD)
    # For example, numbers and arrays will be passed on as into larger numpy arrays
    anlsD = storage[sortd[0]] # A template for anlsD; not necessarily the 0th object
    anlsT = {}
    anlsT["ds_id"] = [None]*numfiles
    for k in anlsD.keys(): # Iterate over key, value pairs for anlsD        
        val = anlsD[k]
        if isinstance(val, numbers.Number): # Single number -> N-element 1D NumPy array
            anlsT[k] = np.zeros(numfiles)
        elif isinstance(val, np.ndarray): # A x B x C x ... Numpy array -> N x A x B x C x ... NumPy array
            anlsT[k] = np.zeros([numfiles] + list(val.shape)) # Consciously choose to keep type as floats, regardless of input type
        else: # Anything else: String, dict, etc. -> N-element list of strings, dicts, etc.
            anlsT[k] = [None]*numfiles

    # Fill the output dictionary anlsT
    for i in range(numfiles):
        anlsD = storage[sortd[i]]
        for k in anlsD:
            anlsT[k][i] = anlsD[k]
    
    if yt.is_root():
        # Dump result to pickled output
        print("All files completed and analysis outputs re-organized. Pickling anlsT...")
        pickle.dump( anlsT, open( os.path.join(outdir, "anlsT.p"), "wb" ) )
        
        # Make the plots
        if callable(plotT):
            print("Making final plots...")
            plotT(anlsT, outdir=outdir)
            
        print("Time-series analysis complete! Outputs stored under " + str(outdir))
    return anlsT


if __name__ == "__main__":
    pass
