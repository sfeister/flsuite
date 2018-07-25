#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
flyt.py: Yt analysis routines for FLASH

Created by Scott Feister on Wed Oct 05 14:51:08 2016
"""

import numpy as np

def get_simdata(ds, lev=None):
    """ Extract a yt covering grid, covering the domain of the entire FLASH simulation
    
    Gives access to underlying data (at fixed resolution) from a FLASH dataset.
    Works in 1D, 2D, or 3D.
        
    Uses yt's "covering grid" method to ensure losslessness.    
    Replaces get_magnetic_field. This will take the AMR Level, and data set,
     and it will return a covering grid in the shape
     of a cube in order to do data analysis.

    Inputs:
        ds:         YT dataset, the output of one timestep e.g. "yt.load(...)"
        lev:        int, level of refinement to access, where 0 is no refinement
                        (yt level 0 matches FLASH lrefine=1)
                        Default value, None, gives maximum available refinement level.
    Outputs:
        cg:         YT covering grid, for the volume of the cube at desired refinement
    
    
    Example usage, 3D Cartesian simulation:
        ds = yt.load('my/folder/mysim_hdf5_plt_cnt_0010')
        # Get a cube out at the yt level 3 of simulation refinement (FLASH lrefine=4)
        cg = get_simdata(ds, lev=3)
        [xmin, ymin, zmin] = cg.left_edge.v  # Lowermost boundaries (out to the cell walls)
        [xmax, ymax, zmax] = cg.right_edge.v  # Uppermost boundaries (out to the cell walls)
        Dens = cg['dens'].v # 3D array of cell density values (axis order: x,y,z)
        
    Example usage, 2D Cartesian simulation:       
        ds = yt.load('my/folder/mysim_hdf5_plt_cnt_0010')
        # Get a rectangle out at the highest possible level of simulation refinement
        cg = get_simdata(ds)
        [xmin, ymin] = cg.left_edge[:2].v
        [xmax, ymax] = cg.right_edge[:2].v
        Dens = np.squeeze(cg['dens'].v) # 2D array of cell density values (axis order: x,y)
        
        fig, ax = plt.subplots()
        ax.pcolorfast([xmin, xmax], [ymin, ymax], Dens.T)
        
    Example usage, 2D Cylindrical simulation:
        ds = yt.load('my/folder/mysim_hdf5_plt_cnt_0010')
        # Get a rectangle out at the highest possible level of simulation refinement
        cg = get_simdata(ds)
        [rmin, zmin] = cg.left_edge[:2].v
        [rmax, zmax] = cg.right_edge[:2].v
        Dens = np.squeeze(cg['dens'].v) # 2D array of cell density values  (axis order: r,z)
        
        fig, ax = plt.subplots()
        ax.pcolorfast([rmin, rmax], [zmin, zmax], Dens.T)
        
    Based on example at http://yt-project.org/doc/examining/low_level_inspection.html#examining-grid-data-in-a-fixed-resolution-array
    Documented, with examples, by J.T. Laune and Scott Feister. Updated July 2018.
    """
    
    if lev is None:
        lev = ds.max_level
    
    if ds.dimensionality == 1:        
        dims = ds.domain_dimensions * 2**np.array([lev, 0, 0])
    elif ds.dimensionality == 2:
        dims = ds.domain_dimensions * 2**np.array([lev, lev, 0])
    else:
        dims = ds.domain_dimensions * 2**lev
            
    cg = ds.covering_grid(level=lev, left_edge=ds.domain_left_edge, dims=dims)
    
    return cg

def get_cube(ds, edge_um=1000, lev=0, center_yt=None):
    """
    Extract a cube of raw data from 3D FLASH yt dataset.
        
    Uses yt's "covering grid" method to ensure losslessness.    
    Replaces get_magnetic_field. This will take the center, edge length,
     AMR Level, and data set, and it will return a covering grid in the shape
     of a cube in order to do data analysis.

    Inputs:
        ds:         YT dataset
        edge_um:    float, cube edge length in microns
        lev:        int, level of refinement to access, where 0 is no refinement (lrefine=1)
        center_yt:  YT array, center position of the cube. If None, domain center is used
    Outputs:
        cg:         YT covering grid, for the volume of the cube at desired refinement
    
    Example usage:
        ds = yt.load('my/folder/mydataset0010')
        # Define cube of density values, offset from domain center in X by 5 microns
        cg = get_cube(ds, edge_um=500, lev=4, center_yt=ds.domain_center + ds.arr([5,0,0],'um')
        X = cg['x'] # 3D array of cell X values
        Y = cg['y'] # 3D array of cell Y values
        Z = cg['z'] # 3D array of cell Z values
        Dens = cg['dens'] # 3D array of cell density values
    
    Documented, with examples, by J.T. Laune and Scott Feister. Updated July 2018.
    """
    
    # Clean up/ check inputs
    if center_yt is None: # Set cube center to domain center, if not provided
        center_yt = ds.domain_center
    elif not isinstance(center_yt, yt.units.yt_array.YTArray): # Verify input center_yt is a YT Array, e.g. not just a NumPy array
        raise TypeError("Input 'center_yt' must be a YT Array (of type yt.units.yt_array.YTArray). E.g. try center_yt = ds.arr([0,0,3],'um')")
     
    # Extract covering-grid cube with given center and edge length
    le = center_yt - ds.arr([edge_um/2, edge_um/2, edge_um/2], 'um') # YT left edge for cube
    re = center_yt + ds.arr([edge_um/2, edge_um/2, edge_um/2], 'um') # YT right edge for cube

    cellwid0 = ds.domain_width/ds.domain_dimensions # Cell width in YT units, at level=0 refinement (lrefine=1)
    
    dims = np.array((re - le)/cellwid0) * 2**lev + 1 # Desired array dimensions ([int, int, int]) of covering grid depends on cell resolution
    cg = ds.covering_grid(level=lev, left_edge=le, dims=dims) # Covering grid cube
    
    return cg

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
