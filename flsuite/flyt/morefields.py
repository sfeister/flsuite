#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
morefields.py: Add a few useful derived fields in future-loaded datasets

Created by Scott Feister around Oct 1 2016
Based on John Zuhone's work / Milad's work / Roman's work in the FLASH YT codebase

Import this module (before loading any data) to add a few derived fields to the yt session.
Additional fields will only be seen in those datasets loaded (e.g. via "yt.load()") AFTER this module has been imported.
Datasets loaded before this module is imported are unaffected.

CHANGELOG:
2016-10-05 Changed name from scottderived.py to morefields.py, moved into the flsuite.flyt sub-module
2016-10-20 Added several more fields, including hooks for kinetic and magnetic energy densities (B**2)

TODO:
* Add magnetic field energy
* Look at what Anthony already wrote; how much of this is already finished there?
"""

from yt import derived_field
#from yt.units import dimensions

@derived_field(name="magtot", units="code_magnetic")
def _magtot(field, data): # Magnitude of the magnetic field vector, returned in code units
    return (data['magx']**2 + data['magy']**2 + data['magz']**2)**0.5

@derived_field(name="magsqr", units="code_magnetic**2")
def _magsqr(field, data): # Square of the magnetic field vector, returned in code units
    return data['magx']**2 + data['magy']**2 + data['magz']**2

@derived_field(name="veltot", units="code_length/code_time")
def _veltot(field, data): # Magnitude of the velocity vector, returned in code units
    return (data['velx']**2 + data['vely']**2 + data['velz']**2)**0.5
    
@derived_field(name="KEdens", units="code_mass/(code_length*code_time**2)")
def _KEdens(field, data): # Classical kinetic energy per volume (KE = 1/2 m v^2), returned in code units
    return 0.5 * data['dens'] * (data['velx']**2 + data['vely']**2 + data['velz']**2)
    
@derived_field(name="cell_KE", units="code_mass*code_length**2/code_time**2")
def _cellKE(field, data): # Classical kinetic energy of each cell (KE = 1/2 m v^2), returned in code units
    return 0.5 * data['cell_mass'] * (data['velx']**2 + data['vely']**2 + data['velz']**2)