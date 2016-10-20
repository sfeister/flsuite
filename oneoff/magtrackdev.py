#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
magtrackdev.py: Scratch space for development of magnetic field tracker

Created by Scott Feister on Wed Oct 19 19:18:17 2016
"""
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    zgv = np.linspace(-0.5,0.5,300)
    velz = np.arcsin(zgv - 0.25) + np.sin(20*zgv)
    
    
    ct = (zgv > -0.1) & (zgv < 0.1)
    zero_crossings = np.where(np.diff(np.signbit(velz)))[0] # Copied from: http://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
    z_crossings = zgv[zero_crossings]
    
    mask_cross = np.append(np.diff(np.signbit(velz)), False) # Boolean array, same size as input, of zero crossings
    if np.sum(mask_cross & ct) > 0:    
        myz = zgv[mask_cross & ct][-1]
    else:
        myz = np.max(zgv[ct])
    
    
    fig = plt.figure(1)
    fig.clear()
    ax = fig.add_subplot(111)
    ax.plot(zgv[ct], velz[ct])
    ax.set_xlabel("Z")
    ax.set_ylabel("Velz")
    ax.axvline(myz, linestyle='--', color='black')
    ax.axhline(0, linestyle='-', color='red')
