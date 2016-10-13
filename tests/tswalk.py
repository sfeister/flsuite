#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
tswalk.py: Throw-away code to try to sort out the 

Created by Scott Feister on Thu Oct 13 15:01:17 2016
"""

import numpy as np
import numbers

def numbers_to_strings(argument):
    # Copied from https://www.pydanny.com/why-doesnt-python-have-switch-case.html
    switcher = {
        0: "zero",
        1: "one",
        2: "two",
    }
    return switcher.get(argument, "nothing")

def allocANLST(anlsD, numfiles):
    """ Intelligently pre-allocate a time-series analysis dictionary using keys
    and values a single-dataset analysis dictionary
    
    anlsD: Single-dataset analysis dictionary (template)
    numfiles: Number of files in the time-series analysis
    anlsT: Time-series dataset analysis dictionary, with keys matching anlsD
    
    For example, numbers and arrays will be passed on as into larger numpy arrays
    """
    
    anlsT = {}
    for k in anlsD.keys(): # Iterate over key, value pairs for anlsD        
        val = anlsD[k]
        if isinstance(val, numbers.Number): # Single number -> N-element 1D NumPy array
            anlsT[k] = np.zeros(numfiles)
        elif isinstance(val, np.ndarray): # A x B x C x ... Numpy array -> N x A x B x C x ... NumPy array
            anlsT[k] = np.zeros([numfiles] + list(val.shape)) # Consciously choose to keep type as floats, regardless of input type
        else: # Anything else: String, dict, etc. -> N-element list of strings, dicts, etc.
            anlsT[k] = [None]*numfiles
    return anlsT
    
if __name__ == "__main__":
    numfiles = 10
    sto = [None]*numfiles
    for i in range(numfiles):
        anlsD = {}
        anlsD["hel"] = "Hello" + str(i)
        anlsD["myarr"] = np.array([i, i + 10, i + 20, i + 30])
        anlsD["myval"] = 20 + i
        anlsD["myfloat"] = 22.0 + i
        anlsD["mynpfloat"] = np.float(22.0)
        sto[i] = anlsD

    # Allocate based on first entry
    anlsD = sto[0]
    anlsT = {}
    mykeys = anlsD.keys()
    for k in mykeys:
        print(k + ": " + str(type(anlsD[k])))
        
        val = anlsD[k]
        if isinstance(val, numbers.Number):
            # Single number
            anlsT[k] = np.zeros(numfiles)
        elif isinstance(val, np.ndarray):
            anlsT[k] = np.zeros([numfiles] + list(val.shape)) # Consciously choose to keep type as floats, regardless of input type
        else:
            #elif isinstance(val, str):
            # String, dict, etc
            anlsT[k] = [None]*numfiles
    
    # Now, fill in the arrays
    for i in range(numfiles):
        anlsD = sto[i]
        for k in mykeys:
            val = anlsD[k]
            anlsT[k][i] = val
