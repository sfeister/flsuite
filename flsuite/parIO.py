#!/bin/python3
# -*- coding: utf-8 -*-
"""
parIO.py: Do input and output on flash.par configuration files

Created by Scott Feister on Fri Sep 16 15:52:46 2016

This is one of the FLASH i/o functions, which are dedicated to parsing and generating inputs and outputs of FLASH.

Changelog:
2016-09-19 110 PM first version of file. For help with Mira re-submission automation.
2016-10-05 Changed name to parIO.py (from parparse1.py) and added to flsuite module
2016-12-15 Updated for single-function restart call, adding function "prepRest()"

Future work:
* Extract log_file and basenm from the flash.par
* Use docopt for filename, log_file etc. input
* Make python-based flash.par reader/writer (general and handling types correctly)

"""

import os
import re
import warnings
import shutil
from datetime import datetime

def getLogNums(basenm, log_file, simdir = '.'):
    """ Get zeroth checkpoint number and first plot output number for a restart
    INPUTS:
    basenm: string, base simulation name specified in flash .par file ('basenm ='), e.g. basenm = "tdyno2016_"
    log_file: string, name of log file with extension as specified in flash .par file ('log_file = '), e.g. log_file = "tdyno2016.log"
    simdir: string, path to simulation output directory, e.g. simdir = "\home\feister\myrun\". Trailing slash optional and has no effect.

    OUTPUTS:
    chknum: integer, checkpoint number to use for restart
    pltnum: integer, first plot number to write during restart run
    
    Example usage:
    chknum, pltnum = getLogNums("tdyno2016_", "tdyno2016.log", "\home\feister\myrun")
    print chknum, pltnum
    """
    
    # Read the log file contents
    with open(os.path.join(simdir, log_file), "r") as f:
        filetxt = f.read() # Read in the whole file, all at once
    # \[IO_write\W*?\] 
        
    matches = re.findall("close: type=[a-z]*? name=" + re.escape(basenm) + "hdf5_([pltchk]*?)(?:_cnt)?_" + "([0-9]+)", filetxt) # Find all instances of closing plotfiles or checkpoints in the log file
    # "matches" is a list of form [('plt','0001'), ...]
    
    # Do a reverse search for plot, then checkpoint
    matches.reverse() # Reverse the matches list
    chknum = None
    pltnum = None
    for m in matches:
        if (not chknum) and m[0] == "chk":
            chknum = int(m[1])
        if chknum and m[0] == "plt":
            pltnum = int(m[1]) + 1
            break
    
    if not chknum:
        raise ValueError('In log file (' + log_file + '), search for final closed checkpoint number failed. Were any checkpoints written?')
    if not pltnum:
        warnings.warn('In log file (' + log_file + '), search for plot file number closed prior to checkpoint ' + str(chknum) + ' failed. Were any plots written? Defaulting to start plot number at 0.')
        pltnum = 0
    if chknum < 0 or pltnum < 0:
        raise ValueError('In log file (' + log_file + '), search ended on a checkpoint number or plot number less than zero. Unsafe to restart!')

    return chknum, pltnum

def writeRePar(chknum, pltnum, simdir = '.', par_file = 'flash.par'):
    """ Modify the (original or restart) flash .par file to reflect a restart from a given checkpoint number, overwriting starting from a given plot number"""

    ## Back up the flash .par (to ".par.bak_DATE") in case we royally screw it up
    shutil.copy(os.path.join(simdir, par_file), os.path.join(simdir, par_file + '.bak_' + "{:%Y%m%d_%H%M}".format(datetime.now())))
    
    ## Read in the .par file contents
    with open(os.path.join(simdir, par_file), "r") as f:
        filetxt = f.read() # Read in the whole file, all at once

    ## Modify the .par contents
    # Step 1: Make sure restart is set to true (restart = .true.)
    filetxt = re.sub("(restart(?:\s*?)=(?:\s*?))(\.[truefalse]+\.)((?:\s*?))", r"\1.true. \3", filetxt)
    # Step 2: Update checkpoint number
    filetxt = re.sub("(checkpointFileNumber(?:\s*?)=(?:\s*?))([0-9]+)((?:\s*?))", r"\g<1>" + str(chknum) + r"\3", filetxt)
    # Step 3: Update plot number
    filetxt = re.sub("(plotFileNumber(?:\s*?)=(?:\s*?))([0-9]+)((?:\s*?))", r"\g<1>" + str(pltnum) + r"\3", filetxt) # The \g<1> is to avoid \1234 regex mistake. http://stackoverflow.com/questions/5984633/python-re-sub-group-number-after-number
    
    ## Write out the modified .par contents (overwrite flash.par)
    with open(os.path.join(simdir, par_file), "w") as f:
        f.write(filetxt)

    return 0

def prepRest(simdir, basenm, lognm, parnm = 'flash.par'):
    """ Prepare a simulation directory to be re-run (restarted)
    High-level function to be called on a simulation directory.
    Details will be replaced down the road.
    
    Example Usage:
        prepRest("\home\feister\myrun", "tdyno2016_", "tdyno2016.log")
    """
    chknum, pltnum = getLogNums(basenm, lognm, simdir) # Analyze the log file for appropriate numbers
    writeRePar(chknum, pltnum, simdir = simdir, par_file = parnm) # Update the flash.par file
    
    return 0
    
def parMod(pardict, simdir='.', par_file='flash.par', outfn=None, bak=True):
    """General modification of a flash.par.
    
    pardict:    Parameters for modification; a python dictionary. A dictionary with variable names as keys and numbers, strings, etc. as values. See following example usage.
    simdir:     Directory of the already-existing .par file
    par_file:   Name (e.g. "flash.par") of the already-existing .par file
    outfn:      Full filepath (e.g. "/home/scott/myrun/flash.par") of the output .par file. If None, overwrite in place.
    bak:        If True, creates a timestamped backup of the input par_file before modification, in its original directory.
    
    Example usage:
    d = {}
    d["basenm"] = '"tdyno2016PI_"' # No spaces allowed here
    d["pi_beamProtonEnergy_1"] = 14.7
    d["pi_beamNumberOfProtons_1"] = 14050821
    d["pi_detectorSideLength_1"] = 9.6 # Detector side length, in cm
    d["pi_beamApertureAngle_1"] = 20.0 # Aperture angle, in degrees
    d["pi_detectorSideTiltingAngle_1"] = 0.0
    d["checkpointFileNumber"] = 220
    d["pi_detectorDist2BeamCapsule_1"] = 40.77 # Distance from beam capsule to detector, in cm
    (d["pi_beamCapsuleX_1"], d["pi_beamCapsuleY_1"], d["pi_beamCapsuleZ_1"]) = (0.0, -1.77, 0.4) # Capsule X,Y,Z position, in cm

    parMod(d, simdir='.', par_file='flash.par') # Updates ./flash.par with the new values
    """
    
    if outfn is None:
        outfn = os.path.join(simdir, par_file) # If no output filename is specified, set it to overwrite the input filename

        ## Back up the flash .par (to ".par.bak_DATE") in case we royally screw it up
    if bak: # By default, input "bak" is True; set "bak=False" to avoid this backup
        shutil.copy(os.path.join(simdir, par_file), os.path.join(simdir, par_file + '.bak_' + "{:%Y%m%d_%H%M}".format(datetime.now())))
    
    ## Read in the .par file contents
    with open(os.path.join(simdir, par_file), "r") as f:
        filetxt = f.read() # Read in the whole file, all at once
    
    for key in pardict:
        # Key will always be a string, of course. TODO: Enforce this somehow else
        filetxt = re.sub("(" + key + "(?:\s*?)=(?:\s*?))(\S+)((?:\s*?))", r"\g<1>" + str(pardict[key]) + r"\3", filetxt) # The \g<1> is to avoid \1234 regex mistake. http://stackoverflow.com/questions/5984633/python-re-sub-group-number-after-number
    
    ## Write out the modified .par contents
    with open(outfn, "w") as f:
        f.write(filetxt)
    
    return 0
 
if __name__ == "__main__":
    pass