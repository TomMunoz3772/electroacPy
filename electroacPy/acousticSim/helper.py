"""
Set of useful functions for manipulating BEM objects and observations

"""

import numpy as np

def sumPressureArray(bemObj, radiatingSurface):
    p_mesh = bemObj.p_mesh
    radSurf_system = bemObj.radiatingElement
    if isinstance(radiatingSurface, str):
        if radiatingSurface == 'all':
            radiatingSurface = radSurf_system

    pressureCoeff = np.zeros([len(bemObj.frequency), 
                              bemObj.spaceP.grid_dof_count], dtype='complex')
    for f in range(len(bemObj.frequency)):
        for j in range(len(radiatingSurface)):
            ind_surface = np.argwhere(radSurf_system == radiatingSurface[j])[0][0]
            pressureCoeff[f, :] += p_mesh[f][ind_surface].coefficients
    return pressureCoeff
