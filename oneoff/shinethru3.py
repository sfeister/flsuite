#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
shinethru3.py: Compute the shine-thru of the TDYNO simulation

Created by Scott Feister around Oct 1 2016

CHANGELOG:
2016-10-05 Copied contents of shinethru2.py into this document, renamed shinethru3.py
"""

## shinethru2.py gets a little more lean in its analysis

print "Importing yt..."
import yt
#from yt import YTArray
print "Performing other imports..."
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


print "Beginning script."

datdir = r'/gpfs/mira-fs0/projects/CosmicLaser/tzeferac/NIF/TDYNO_BAND2'

for i in range(50, 272, 272): # Iterate over plot files (time steps)
	# Build the HDF5 filename
	name = r'tdyno2016_hdf5_plt_cnt_' + str(i).zfill(4)
	print "WORKING ON:", name
	fn = os.path.join(datdir, name)
	
	# Load the dataset
	ds = yt.load(fn) # Load the flash plot file

	## Make a lineout plot
	# x lineout; cut through y = z = 0
	ray = ds.ortho_ray(2, (ds.domain_center[0], ds.domain_center[1]))
	srt = np.argsort(ray['z'])
	
	# Plot lineout in matplotlib
	plt.figure(1)
	plt.clf()
	plt.subplot(211)
	plt.plot(np.array(ray['z'][srt]), np.array(ray['density'][srt]))
	plt.ylabel('density')
	plt.subplot(212)
	plt.plot(np.array(ray['z'][srt]), np.array(ray['depo'][srt]))
	plt.xlabel('z')
	plt.ylabel('deposition')

	plt.savefig("den_depo_xsweep.png")
	
	## Make a slightly well logically-defined frb
	# User specifications
	xyzlims = ds.arr([[-2, 2], [-2, 2], [-2, 10]], 'mm') #x, y, z limits. [[xmin, xmax], [ymin, ymax], [zmin, zmax]]. Also defines units for plot.
	axis = 0 # Slicing axis
	pixmax = 1000 # Number of pixels along the longest axis
	fldnm = 'depo' # Field name (Could this be a list? Not sure!)
	
	# Extract a relevant sub-region, and assign some basic values
	centr = np.mean(xyzlims,1) # Center for each dimension. Equivalent to centr = myreg.center()
	widt = np.ptp(xyzlims,1) # Width for each dimension. Equivalent to widt = myreg.left_edge + (myreg.right_edge - myreg.left_edge)/2
	myreg = ds.region(centr, xyzlims[:,0], xyzlims[:,1], fields=(fldnm)) # Sub-region of data, only in the interesting parts.
	
	# Perform a 2D projection on the region
	proj = myreg.integrate(fldnm, axis) # Do the projection (why do we even bother with a field name here?)
	
	# Convert 2D projection to 2D fixed-resolution-buffer (frb)
	wid2Dtmp1 = np.delete(widt.v, axis) # Width for each ~plotting~ dimension (projection dimension removed); Two-element NDArray
	res2D = ((wid2Dtmp1 / np.max(wid2Dtmp1)) * pixmax)[::-1] # Resolution for each plotting dimension. NDArray. Specify the biggest axis as pixmax pixels, and the other one appropriately scaled. Resolution two-element list has to be reversed (the [::-1]) relative to width two-elements. Two-element NDArray
	wid2Dtmp2 = np.array([np.delete(widt.v, axis), [widt.units]*2]).T # Widths of plotting x,y axes, formatted as numpy array with units; 2x2 NDArray e.g. [[11.0, 'mm'], [12.0, 'mm']]
	wid2D = tuple(map(tuple, wid2Dtmp2)) # Widths of plotting x,y axes, formatted as tuple of tuples; e.g. ((11.0, 'mm'), (12.0, 'mm'))
	del wid2Dtmp1, wid2Dtmp2

	frb = proj.to_frb(wid2D, res2D, center=centr) # Cast to fixed-resolution buffer (this is actually usable with any field name)

	# Make a 2D colormap plot using frb
	lims2D = np.delete(xyzlims, axis, axis=0) # Limits for plotting x,y axes. Keeps as a YTArray, but now it is 2x2 like YTArray [[-2,2], [-10,10]]
	lbl2D = np.delete(np.array(["X","Y","Z"]), axis) # Get the labels of the slice by removing the slicing axis

	[[xmin, xmax], [ymin, ymax]] = lims2D.v
	units = str(lims2D.units)
	
	plt.figure(4)
	plt.clf()
	plt.imshow(
		np.log10(np.array(frb[fldnm])), 
		cmap='viridis', 
		extent=[xmin, xmax, ymin, ymax],
		)
	#plt.title("Log10 Projection of " + fldnm + " (" + str(frb[fldnm].units) + ")")
	plt.title("Log10 Integ(" + fldnm + ")")
	plt.xlabel(lbl2D[0] + " (" + units + ")")
	plt.ylabel(lbl2D[1] + " (" + units + ")")
	plt.colorbar()
	plt.savefig('my_fav_figure.png')
	
	## Try a simple cylinder projection
	centr = ds.domain_center
	rad = ds.quan(1, 'mm')
	hgt = ds.quan(10, 'mm')
	mycyl = ds.disk(centr, [0, 0, 1], rad, hgt)
	#proj2 = mycyl.integrate(fldnm, axis)
	p = yt.ProjectionPlot(ds, 'x', "density", method='integrate', width = [(10,'mm'),(1,'mm')], data_source=mycyl)
	p.set_log('density', False)
	p.save()
	
	## Try a cylinder line plot
	ds.add_field(("misc", "ones"), function=lambda x, y: 1.0) # Add this derived field, which is always dimensionless one, for normalization
	centr = ds.domain_center
	#centr = ds.arr([0.0, 0.0, 0.0], 'mm')
	rad = ds.quan(0.05, 'mm')
	hgt = ds.quan(6, 'mm') # Note: This 'hgt' is actually HALF the cylinder's total height, per yt spec
	mycyl = ds.disk(centr, [0, 0, 1], rad, hgt)
	projA = mycyl.integrate('dens', 'x')
	projB = mycyl.integrate('dens', 'y')

	# Convert 2D projection to 2D fixed-resolution-buffer (frb)
	res2D = np.array([2000, 1000]) # [zaxis pixels, xaxis pixels]
	wid2D = ((2*rad.v, str(rad.units)), (2*hgt.v, str(hgt.units))) # e.g. ((11.0, 'mm'), (12.0, 'mm'))
	frbA = projA.to_frb(wid2D, res2D, center=centr) # Cast to fixed-resolution buffer
	frbB = projB.to_frb(wid2D[::-1], res2D[::-1], center=centr) # Cast to fixed-resolution buffer


	plt.figure(1)
	plt.clf()
	plt.imshow(frbA["dens"].v)
	plt.savefig('dbgA.png')
	plt.clf()
	plt.imshow(frbB["dens"].v)
	plt.savefig('dbgB.png')


	lineout = (np.nanmean(frbA["dens"]/frbA["ones"], 1) + np.nanmean(frbB["dens"]/frbB["ones"], 0))/2
	zgv = np.linspace(centr[2] - hgt, centr[2] + hgt, res2D[0])


	plt.figure(5)
	plt.clf()
	plt.plot(zgv.in_units('mm'), lineout)
	plt.scatter(ray['z'][srt].in_units('mm').v, ray['density'][srt].v)
	plt.xlabel("Z (mm)")
	plt.ylabel("Line-width-averaged density (" + str(lineout.units) + ")")
	plt.xlim(0.2, 0.7)
	plt.ylim(-0.1, 0.7)
	plt.savefig('lineplot_dens.png')

	plt.figure(6)
	plt.clf()
	plt.plot(ray['z'][srt].in_units('mm').v, ray['velz'][srt].v)
	plt.xlabel("Z (mm)")
	plt.ylabel("Fluid velocity")
	plt.xlim(-2, 2)
	plt.savefig('lineplot_velz.png')
