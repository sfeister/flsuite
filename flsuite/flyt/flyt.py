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
    
if __name__ == "__main__":
    pass
