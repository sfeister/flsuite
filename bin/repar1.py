#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
repar1.py: Reformat flash.par file for a quick and dirty restart, using log file to set checkpoint/plot numbers

Created by Scott Feister on Thu Dec 15 15:19:34 2016

EXAMPLE USE CASE:
    cd ~/my/sim/run/dir
    mpirun flash4
    [Runs simulation... then stops]
    python ~/my/python/dir/repar1.py
    mpirun flash4
    
CHANGELOG:
2016-12-15  Created file, based on extracting elements from parIO.py (in Scott's flsuite module) -SF
            Added an attempt to get the basenm and lognm from the flash.par file, for better portability -SF

TODO:
* Check that flash.par, lognm, and basenm files exist (throw errors if not)

"""

import os
import re
import warnings
import shutil
from datetime import datetime

if __name__ == "__main__":
    ### USER: DEFINE YOUR PARAMETERS HERE
    simdir = '.' # Assume the working directory contains simulation outputs
    parnm = 'flash.par' # Assume the .par is called "flash.par"

    ## Read in the .par file contents (to get basename and construct output log filename)
    with open(os.path.join(simdir, parnm), "r") as f:
        partxt = f.read() # Read in the whole file, all at once
    basenm = re.findall("basenm(?:\s*?)=(?:\s*?)\"(.*)\"", partxt)[0] # Extract basename from flash.par. Could probably do this better, but hey, it works (for now!)

    if basenm.endswith("_"):
        basenm_shrt = basenm[:-1] # No underscore, e.g. "tdyno2016"
        basenm_lng = basenm # Has underscore, e.g. "tdyno2016_"
    else:
        basenm_shrt = basenm # No underscore, e.g. "tdyno2016"
        basenm_lng = basenm + '_' # Has underscore, e.g. "tdyno2016_"

    lognm = basenm_shrt + '.log' # Guesses at .log filename, e.g. "tdyno2016.log"
    
    # Read the log file contents
    with open(os.path.join(simdir, lognm), "r") as f:
        logtxt = f.read() # Read in the whole file, all at once
        
    matches = re.findall("close: type=[a-z]*? name=" + re.escape(basenm_lng) + "hdf5_([pltchk]*?)(?:_cnt)?_" + "([0-9]+)", logtxt) # Find all instances of closing plotfiles or checkpoints in the log file
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
        raise ValueError('In log file (' + lognm + '), search for final closed checkpoint number failed. Were any checkpoints written?')
    if not pltnum:
        warnings.warn('In log file (' + lognm + '), search for plot file number closed prior to checkpoint ' + str(chknum) + ' failed. Were any plots written? Defaulting to start plot number at 0.')
        pltnum = 0
    if chknum < 0 or pltnum < 0:
        raise ValueError('In log file (' + lognm + '), search ended on a checkpoint number or plot number less than zero. Unsafe to restart!')

    ### Modify the (original or restart) flash .par file to reflect a restart from a given checkpoint number, overwriting starting from a given plot number

    ## Back up the flash .par (to ".par.bak_DATE") in case we royally screw it up
    shutil.copy(os.path.join(simdir, parnm), os.path.join(simdir, parnm + '.bak_' + "{:%Y%m%d_%H%M}".format(datetime.now())))
    
    ## Modify the .par contents in RAM
    # Step 1: Make sure restart is set to true (restart = .true.)
    partxt = re.sub("(restart(?:\s*?)=(?:\s*?))(\.[truefalse]+\.)((?:\s*?))", r"\1.true. \3", partxt)
    # Step 2: Update checkpoint number
    partxt = re.sub("(checkpointFileNumber(?:\s*?)=(?:\s*?))([0-9]+)((?:\s*?))", r"\g<1>" + str(chknum) + r"\3", partxt)
    # Step 3: Update plot number
    partxt = re.sub("(plotFileNumber(?:\s*?)=(?:\s*?))([0-9]+)((?:\s*?))", r"\g<1>" + str(pltnum) + r"\3", partxt) # Note: the \g<1> is to avoid \1234 regex mistake. http://stackoverflow.com/questions/5984633/python-re-sub-group-number-after-number
    
    ## Write out the modified .par contents (overwrite flash.par)
    with open(os.path.join(simdir, parnm), "w") as f:
        f.write(partxt)
    
    print("Updated flash.par to restart from checkpoint #" + str(chknum).zfill(4) + " and plot #" + str(pltnum).zfill(4) + ".")