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
2016-11-18 Added abar (average atomic mass), zbar (average ionization state) based on Petros' definitions of 'ye  ', 'sumy' -ScottFeister
           Added sound speed, mach number
2016-11-29 Added nion, nele from FLASH user manual. And Reynolds number, etc.; wait... were these already given by derived fields of the yt modules?
           
TODO:
* Test validity of sound speed, Mach number derived fields
* Fix temperature units (they are in Kelvin) -- enable conversion to eV
* Add nion from Anthony's equation: return data['dens'] * data['sumy'] * 6.0221409E23
"""

from yt import derived_field
from yt.units import centimeter, second, gram
import numpy as np

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
    
@derived_field(name="abar", units="dimensionless")
def _abar(field, data): # Average atomic mass (amu) of an ion
    return 1.0 / data['sumy']
    
@derived_field(name="zbar", units="dimensionless")
def _zbar(field, data): # Average ionization of an ion
    return data['ye  '] / data['sumy']

@derived_field(name="cs", units="cm/s")
def _cs(field, data): # Sound speed, calculated from Petros' 2016 paper table; converted temps to eV
    return (centimeter / second) * 9.80e5 * (data['zbar'] * data['tele'].to_equivalent('eV','thermal').v + (5./3.) * data['tion'].v/11604.525)**0.5 / data['abar']**0.5

@derived_field(name="coulog", units="dimensionless")
def _coulog(field, data): # Coulomb logarithm, ln(Gamma). Includes Kelvin conversion to eV 11604.525
    # 23.5 - np.log(data['nele'].in_units('1/cm**3').v**0.5 * (data['tele'].v/11604.525)**(-1.25)) - np.sqrt(1.0e-5 + (np.log(data['tele'].v/11604.525) - 2)**2/16.0)
    dimensionless = (gram/gram) # Silly hack needed for yt; need to multiply by dimensionless or gets a units error; not sure what the deal is -ScottFeister
    return (23.5 - np.log(data['nele'].in_units('1/cm**3').v**0.5 * (data['tele'].to_equivalent('eV','thermal').v)**(-1.25)) - np.sqrt(1.0e-5 + (np.log(data['tele'].to_equivalent('eV','thermal').v) - 2)**2/16.0)) * dimensionless

@derived_field(name="mach", units="dimensionless")
def _mach(field, data): # Mach number (positive-definite)
    return data['veltot'] / data['cs']

@derived_field(name="nion", units="1/code_length**3")
def _nion(field, data): # Ion number density; see "Examining LaserSlab Simulation Output" gray box in FLASH User Manual
    return data['dens'] * data['sumy'] * (6.0221409E23/gram) # (1 g = 6.022E23 amu)

@derived_field(name="nele", units="1/code_length**3")
def _nele(field, data): # Electron number density; see "Examining LaserSlab Simulation Output" gray box in FLASH User Manual
    return data['zbar'] * data['nion']

@derived_field(name="iimfp", units="cm")
def _iimfp(field, data): # Ion-ion mean-free-path, in cm
    return 2.88e13 * centimeter * (data['tion'].to_equivalent('eV','thermal').v)**2 / (data['zbar']**4 * data['nion'].in_units('1/cm**3').v * data['coulog'])

@derived_field(name="RmperL", units="1/cm")
def _RmperL(field, data): # Magnetic Reynolds number per length scale L; assuming Rm = uL/eta, and RmperL = u/eta
    eta = (3.2e5 * centimeter**2 / second) * data['zbar'] * data['coulog'] / (data['tele'].to_equivalent('eV','thermal').v)**1.5
    return data['veltot'] / eta

@derived_field(name="ReperL", units="1/cm")
def _ReperL(field, data): # Reynolds number per length scale L; assuming Re = uL/nu, and ReperL = u/nu
    nu = (1.92e19 * centimeter**2 / second) * (data['tion'].to_equivalent('eV','thermal').v)**2.5 / (data['abar']**0.5 * data['zbar']**4 * data['nion'].in_units('1/cm**3').v * data['coulog'])
    return data['veltot'] / nu

@derived_field(name="Pm", units="dimensionless")
def _Pm(field, data): # Prandtl number, for simple case of fixed length as in Petros' paper Table II. Pm = Rm/Re; shared length scales
    return data['RmperL']/data['ReperL']

