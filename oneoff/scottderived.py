#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
scottderived.py: Import this module to add a few magnetic-field-based derived fields to the yt session

Created by Scott Feister around Oct 1 2016

"""

# Add a few derived fields for flash.

from yt import derived_field
from yt.units import dimensions

@derived_field(name="magtot", units="code_magnetic")
def _magtot(field, data): # Magnitude of the magnetic field vector, returned in code units (nah... having issues here.)
    return (data['magx']**2 + data['magy']**2 + data['magz']**2)**0.5
	
@derived_field(name="veltot", units="auto", dimensions=dimensions.velocity)
def _veltot(field, data): # Magnitude of the velocity vector, returned in code units
    return (data['velx']**2 + data['vely']**2 + data['velz']**2)**0.5