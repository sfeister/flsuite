#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
flyt.py: Yt analysis routines for FLASH

Created by Scott Feister on Wed Oct 05 14:51:08 2016
"""

import numpy as np

# Would like to get this to taking multiple fields, methods, dimensions
def lineout(ds, field='density', method='ray'):
    """Return a lineout of a yt dataset"""
    if method == 'ray': # Do a ray-style lineout. Spacing not guaranteed.
        # Ray must be sorted. See note on ray data ordering: http://yt-project.org/doc/faq/index.html#ray-data-ordering
        ray = ds.ortho_ray(2, (ds.domain_center[0], ds.domain_center[1]))
        srt = np.argsort(ray['t'])
    
        zvals = np.array(ray['z'][srt])
        fldvals = np.array(ray[field][srt])
    else:
        return 1
    
    return zvals, fldvals
    
if __name__ == "__main__":
    pass
