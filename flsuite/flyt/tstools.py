#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
tstools.py: Tools for manipulating YT Time Series (bundles of FLASH datasets)

Created by Scott Feister on Wed Oct 05 18:18:14 2016

Example usage:
import flsuite.flyt.tstools as tst

CHANGELOG:
2016-10-05 Created tstools.py, filled with some functions from tstest3.py
2016-10-13 Greatly expanded the scope, automating post-analysis.
2018-03-15 Added support for Python3
2018-07-26 Major reorganization for publication online. -SF

TODO:
* Tie together plotT and anlsD, perhaps making a class!
* Enable further modularity to allow multiple instances of anlsD, plotT to be included (e.g. allows a list of functions for each)
   Possibly, even have the anlsD allow anlsT inputs (for auto-scaling, etc.). Or, allow classes manually setting parameters like mins, maxes, etc. Could become a monstrosity if that, though!

############# EXAMPLE SCRIPT #############

# Note: This script can be run in serial (python myscript.py) or in parallel (mpirun -n 4 python myscript.py)

import os
import yt
import numpy as np
import flsuite.sftools as sf
import flsuite.flyt.tstools as tst

def anlzD(ds, anlsD, outdir):
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

def plotT(anlsT, outdir):
    plt.figure(1)
    plt.clf()
    plt.plot(anlsT['times_ns'], anlsT['Bav_code']*np.sqrt(4*np.pi)*1e-3, label='All space')
    plt.xlabel("Simulation time (ns)")
    plt.ylabel("Mean magnetic field strength (kG)")
    plt.title("TDYNO_BAND, Magnetic field through time")
    plt.legend(loc='upper left')
    plt.savefig(os.path.join(outdir,'shinethru1.png'), dpi=300)

    
ts = yt.load("/my/flash/outs/lasslab_hdf_plt_cnt_????")
anlsT = tst.anlzT(ts, anlzD, outdir="/my/plot/path", plotT=plotT, trypkl=True)

######### END OF EXAMPLE SCRIPT #########
"""

import os
import sys
import numpy as np
if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import pickle
else:
    # Python 2 code in this block
    import cPickle as pickle
import yt
yt.enable_parallelism() # Tap into yt's mpi4py parallelism (e.g. now can call via mpirun -np 10 python <blah>.py)
yt.funcs.mylog.setLevel(30) # This sets the output notification threshold to 30, WARNING. Default is 20, INFO.
from mpi4py import MPI
import numbers
import inspect

comm = MPI.COMM_WORLD
rank = comm.rank

# TODO: Incorporate this into the new analysis schema
def tsinf(ts):
	""" Walks through the time series and returns a dict with some basic info like time steps """
	numfiles = len(ts)

	ts_inf = {}
	ts_inf['times_ns'] = np.zeros(numfiles)
	
	for i in range(numfiles):
		ts_inf['times_ns'][i] = ts[i].current_time.in_units('ns').v

	return ts_inf

# TODO: Better documentation
# TODO: Homogenize string outputs from anlsD such that they become numpy arrays (?)
def anlzT(ts, anlzD, outdir='.', plotT=None, trypkl=True):
    """(*AnaLyZe Time series*) Walks through the time series and does a custom analysis "anlzone(ds)" IN PARALLEL on all datasets in time series
    Follows walkthrough at: http://yt-project.org/doc/analyzing/parallel_computation.html?highlight=parallel#parallelization-over-multiple-datasets-including-time-series
    
    Inputs for this function:
        ts:  Yt time series
        anlzD: Function defined with required inputs (ds, anlsD, outdir). Runs on every dataset. Can be written to modify "anlsD" dict in-place by adding values and plot images into "outdir".
        outdir: A directory into which to save your plots and anlsT.p pickle file.
        plotT: (Optional) Function defined with inputs (anlsT, outdir). Uses anlsT to make plots into directory 'outdir'. Run once at end of analysis.
        trypkl: If True (default), will attempt to load the "anlsT.p" pickle file in the directory if it already exists (skipping HDF5 re-analysis). If False, will ignore and overwrite this file.
    
    More info on dictionaries used of this function:
        anlsD: Dictionary compiling output from a single dataset.
        anlsT: Dictionary compiling outputs from time-series, keys matching those of anlsD.

    With trypkl=False, will always re-analyze the HDF5 files.
    With trypkl=True, will nonetheless check to see if pickle file exists, then unpickle and check if any obvious changes have occurred in the anlzD or HDF5 filename list before continuing onto the plotT function.
    
    If run in parallel, only the root will output the correct copy of anlsT; all others will output {}.
    
    Examples of these functions, which would be included in your code, are included in the help string at the top of this file.
    """
    
    anlsT = {}
    pfile = os.path.join(outdir, "anlsT.p")
    
    usepfile = trypkl and os.path.exists(pfile) # True or False, should we try and use the pickle file anlsT.p? Depends on option 'trypkl' and if that file exists.
    if usepfile: # Attempt to load anlsT from anlsT.p
        if yt.is_root():
            print("Found pickle file in output directory: " + str(pfile))            
            print("Unpickling anlsT.p.")
            
            try:
                anlsT = pickle.load( open( pfile, "rb" ) )                    
                print("Successfully unpickled anlsT.p.")
            except Exception as e:
                print("Caught this error: " + repr(e))
                raise Exception("Unpickling anlsT.p failed. Perhaps it was pickled using another version of Python? Delete file or set anlzT option 'trypkl=False' and run again.")
            
            if not np.array_equal(np.array(anlsT['ts_outputs']), np.array(ts.outputs)): # Are ts.outputs and ts_outputs the same? Don't assume lists or numpy arrays as input, force to numpy arrays for comparisons
                print("HDF5 filename list changed, so forcing a complete re-analysis.")
                usepfile = False
            elif anlsT['anlzD_co_code'] != anlzD.__code__.co_code: # Check if bytecode of anlzD has changed (will catch some but not all changes, e.g. will catch new lines of code but miss changes in strings or numbers)
                print("Bytecode changes detected in anlzD, so forcing a complete re-analysis.")
                usepfile = False
            else:
                try:
                    if anlsT['anlzD_source'] != inspect.getsource(anlzD):
                        print("Source code changes detected in anlzD, so forcing a complete re-analysis.")
                        usepfile = False
                except:
                    pass
                    
    usepfile = comm.bcast(usepfile, root=0)
                
    if usepfile: # Use the pickle file directly
        if yt.is_root():
            print("Skipping anlzD function, analysis of HDF5 files (set anlzT option 'trypkl=False' to force re-analysis).")
    else: # Use anlzD to generate anlsT from the HDF5 files
        numfiles = len(ts)

        if yt.is_root():
            print("Number of files to analyze:" + str(numfiles))
            
        ## Analyze datasets using yt parallelization
        storage = {} # Define an empty storage dictionary for collecting information in parallel through processing

        for sto, ds in ts.piter(storage=storage): # Iterate over datasets ds in time-series ts
            print("Reading file:" + str(ds))
            anlsD = {}
            anlsD["ds_id"] = str(ds)
            anlsD['ds_time'] = ds.current_time.v
            anlzD(ds, anlsD, outdir) # It would be nice not to require the outdir argument
            
            if not isinstance(anlsD, dict): # If for some reason anlzD mangled the anlsD array into something other than a dictionary
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
        
        anlsT['anlzD_co_code'] = anlzD.__code__.co_code # Store bytecode instructions for the function anlzD for later comparison
        anlsT['ts_outputs'] = ts.outputs # Store a list of the input files for later comparison
        try:
            anlsT['anlzD_source'] = inspect.getsource(anlzD) # Store the source code of anlzD for comparison (only works if anlzD is defined in a file)
        except:
            pass

        if yt.is_root():
            # Dump result to pickled output
            print("All files completed and analysis outputs re-organized. Pickling anlsT...")
            pickle.dump( anlsT, open( pfile, "wb" ) )
    

    # Make the plots
    if yt.is_root():
        if callable(plotT):
            print("Making final plots (calling plotT function)...")
            plotT(anlsT, outdir)
            
        print("Time-series analysis complete! Outputs stored under " + str(outdir))
        
    return anlsT


if __name__ == "__main__":
    pass
