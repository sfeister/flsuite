#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
flyt.py: Yt analysis routines for FLASH

Created by Scott Feister on Wed Oct 05 14:51:08 2016
"""

import numpy as np

# Would like to get this to taking multiple fields, methods, dimensions. Like to be able to lineout other than center as well
def lineout(ds, field='density', method='ray', npts=1000, axis=0):
    """Return a lineout of a yt dataset, along z, through X=Y=0
    
    npts: Number of points in the lineout
    """
    if method == 'ray': # Do a ray-style lineout. Spacing not guaranteed.
        # Ray must be sorted. See note on ray data ordering: http://yt-project.org/doc/faq/index.html#ray-data-ordering
        #axstr = 'xyz'[axis]
        ax1 = axis
        ax2 = np.delete([0,1,2],ax1)[0]
        ax3 = np.delete([0,1,2],ax1)[1]
        
        ray = ds.ortho_ray(ax1, (ds.domain_center[ax2], ds.domain_center[ax3]))
        srt = np.argsort(ray['xyz'[ax1]])
    
        rayx = np.array(ray['xyz'[ax1]][srt])
        rayfld = np.array(ray[field][srt])

        xgv = np.linspace(np.min(rayx), np.max(rayx), 200)
        fldgv = np.interp(xgv, rayx, rayfld)
    else:
        return 1
    
    return xgv, fldgv

def lineout2(ds, flds=['density'], method='ray', npts=1000, axis=0):
    """Return a lineout of a yt dataset, along z, through X=Y=0
    
    npts: Number of points in the lineout
    
    Expects multiple fields. Also, outputs YT arrays rather than pure NP arrays.
    """
    if method == 'ray': # Do a ray-style lineout. Spacing not guaranteed.
        # Ray must be sorted. See note on ray data ordering: http://yt-project.org/doc/faq/index.html#ray-data-ordering
        #axstr = 'xyz'[axis]
        ax1 = axis
        ax2 = np.delete([0,1,2],ax1)[0]
        ax3 = np.delete([0,1,2],ax1)[1]
        
        ray = ds.ortho_ray(ax1, (ds.domain_center[ax2], ds.domain_center[ax3]))
        srt = np.argsort(ray['xyz'[ax1]])
    
        rayx = ray['xyz'[ax1]][srt] # YT array
        xgv = np.linspace(rayx.min(), rayx.max(), 200) # linspace retains the units successfully

        fldgv = {}
        for fld in flds:
            rayfld = ray[fld][srt]
            fldgv[fld] = ds.arr(np.interp(xgv, rayx, rayfld), rayfld.units) # xgv and rayx need to be in same units, as they are always. np interp loses the units, so add them back in
    else:
        return 1
    
    return xgv, fldgv

# TODO: Arbitrary axis lineout, not just Z
def fatline(ds, rad_mm, fld='density', npts=1000):
    """Return a lineout of a yt dataset, along z, through X=Y=0 but averaged over X,Y within a radius
    
    Uses a cylinder region rather than ray to achieve averaging
    
    """
    if 'ones' not in [x[1] for x in ds.field_list]:
        ds.add_field(("misc", "ones"), function=lambda x, y: 1.0) # Add this derived field, which is always dimensionless one, for normalization
        
    centr = ds.domain_center
    rad = ds.quan(rad_mm, 'mm') # YT array, in mm units
    hgt = ds.domain_width[2]/2 # YT array, in code_length units. Note: This 'hgt' is actually HALF the cylinder's total height, per yt spec for subsequent call to ds.disk()
    mycyl = ds.disk(centr, [0, 0, 1], rad, hgt)
    projA = mycyl.integrate(fld, 'x')
    projB = mycyl.integrate(fld, 'y')

    # Convert 2D projection to 2D fixed-resolution-buffer (frb)
    res2D = np.array([npts, 1200]) # [zaxis pixels, xaxis pixels]
    wid2D = ((2*rad.v, str(rad.units)), (2*hgt.v, str(hgt.units))) # e.g. ((11.0, 'mm'), (12.0, 'mm'))
    frbA = projA.to_frb(wid2D, res2D, center=centr) # Cast to fixed-resolution buffer
    frbB = projB.to_frb(wid2D[::-1], res2D[::-1], center=centr) # Cast to fixed-resolution buffer

    fldgv = (np.nanmean(frbA[fld]/frbA["ones"], 1) + np.nanmean(frbB[fld]/frbB["ones"], 0))/2
    zgv = np.linspace(centr[2] - hgt, centr[2] + hgt, res2D[0])

    return zgv, fldgv
    
def fatline2(ds, rad_mm, flds=['density'], npts=1000):
    """Return a lineout of a yt dataset, along z, through X=Y=0 but averaged over X,Y within a radius
    
    Uses a cylinder region rather than ray to achieve averaging
    
    Expects multiple fields; returns dict of field grid vectors.
    Also, outputs YT arrays (keep units)
    """
    if 'ones' not in [x[1] for x in ds.field_list]:
        ds.add_field(("misc", "ones"), function=lambda x, y: 1.0) # Add this derived field, which is always dimensionless one, for normalization
        
    centr = ds.domain_center
    rad = ds.quan(rad_mm, 'mm') # YT array, in mm units
    hgt = ds.domain_width[2]/2 # YT array, in code_length units. Note: This 'hgt' is actually HALF the cylinder's total height, per yt spec for subsequent call to ds.disk()
    mycyl = ds.disk(centr, [0, 0, 1], rad, hgt)
    projA = mycyl.integrate(flds[0], 'x')
    projB = mycyl.integrate(flds[0], 'y')

    # Convert 2D projection to 2D fixed-resolution-buffer (frb)
    res2D = np.array([npts, 1200]) # [zaxis pixels, xaxis pixels]
    wid2D = ((2*rad.v, str(rad.units)), (2*hgt.v, str(hgt.units))) # e.g. ((11.0, 'mm'), (12.0, 'mm'))
    frbA = projA.to_frb(wid2D, res2D, center=centr) # Cast to fixed-resolution buffer
    frbB = projB.to_frb(wid2D[::-1], res2D[::-1], center=centr) # Cast to fixed-resolution buffer

    zgv = np.linspace(centr[2] - hgt, centr[2] + hgt, res2D[0])
    fldgv = {}
    for fld in flds:
        fldgv[fld] = (np.nanmean(frbA[fld]/frbA["ones"], 1) + np.nanmean(frbB[fld]/frbB["ones"], 0))/2
        
    return zgv, fldgv

def zsplits(ds, flds=['density'], npts=21):
    """Return a lineout of a yt dataset, along z, through X=Y=0 but averaged over X,Y within a cube edge length
    
    Splits region into a series of cubes and takes their averages
    
    Expects multiple fields; returns dict of field grid vectors.
    Also, outputs YT arrays (keep units)
    """
    # Define the centers
    zmin = ds.domain_left_edge[2]
    zmax = ds.domain_right_edge[2]
    
    dz = (zmax - zmin)/npts
    zgv = np.linspace(zmin + dz, zmax - dz, npts)
    fldgv = {}
    
    for i in range(npts):
        # For each zvalue, define a cube centered on (X, Y, Z) = [0, 0, zcent] and take its average
        cent = ds.arr([ds.domain_center[0].v, ds.domain_center[1].v, zgv[i].v], 'code_length') # Typically, centered at (0, 0, zcent)
        le = ds.arr([ds.domain_left_edge[0].v, ds.domain_left_edge[1].v, (zgv[i] - dz/2.0).v], 'code_length') # Sub-region left edge
        re = ds.arr([ds.domain_right_edge[0].v, ds.domain_right_edge[1].v, (zgv[i] + dz/2.0).v], 'code_length') # Sub-region right edge
        reg = ds.region(cent, le, re) # Sub-region rectangular prism
                
        for fld in flds:
            myval = reg.mean(fld)
            if i < 1: # First iteration, allocate arrays
                fldgv[fld] = ds.arr(np.zeros([npts]), myval.units)
            fldgv[fld][i] = myval
            reg.clear_data() # Free up memory after each field; probably good practice given high-res 3D datasets
        
    return zgv, fldgv


if __name__ == "__main__":
    pass
