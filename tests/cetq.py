#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
cetq.py: Check the Cetus queue using qstat.

Created by Scott Feister on Wed May 31 15:31:50 2017
"""

from subprocess import check_output
import os

def get_jobs(user = None):
    if user is None:
        jobsblob = check_output(["qstat", "-lf"])
    else:
        jobsblob = check_output(["qstat", "-lf", "-u", user])    
    jobslist = jobsblob.split(os.linesep + os.linesep)
    
    jobsdict = {}
    
    for jobstr in jobslist:
        linelist = jobstr.split(os.linesep)
        if len(linelist) > 2: # Make sure this isn't a blank line
            jobdict = {}
            
            for line in linelist:
               parts = line.split(":")
               jobdict[parts[0].strip()] = parts[1].strip()
            jobsdict[jobdict["JobID"]] = jobdict

    return jobsdict
    
if __name__ == "__main__":
    jobsdict = get_jobs(user="sfeister")
    print jobsdict.keys()