"""
Save and Load projects
"""
import numpy as np
import os
from shutil import copy2
from os.path import join
from numpy import asanyarray as array
from electroacPy import loudspeakerSystem
from electroacPy.acousticSim.bem import bem
from electroacPy.acousticSim.pointSource import pointSource, pointSourceBEM
from electroacPy.acousticSim.evaluations import evaluations as evs
import bempp_cl.api

def save(projectPath, system):
    """
    Save system objects in NPZ files.

    Parameters
    ----------
    projectPath : TYPE
        DESCRIPTION.
    loudspeakerSystem : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    if not os.path.exists(projectPath):
        os.mkdir(projectPath)
    
    
    # main variables
    frequency   = system.frequency
    c           = system.c
    rho         = system.rho
    radiator_id = system.radiator_id
    
    # save main variables + LEM objects
    np.savez(join(projectPath, 'LEM'),
             frequency   = frequency,
             driver      = system.driver,
             laser       = system.vibrometry,
             crossover   = system.crossover,
             enclosure   = system.enclosure,
             radiator_id = radiator_id,
             c           = c,
             rho         = rho)
    
    # save BEM parameters
    

    
    
    
    
    return None

def load(pathToProject):
    """
    Load loudspeaker simulation project from directory path.

    Parameters
    ----------
    pathToProject: str,
        path to directory where npz files are stored

    Returns
    -------
    LS: loudspeakerSystem object,
        Study, evaluations and LEM setup.
    """

    return None


#%% HELPERS
def storePressureMeshResults_PSADM(acoustic_study):
    study = acoustic_study
    Nfft = len(study.frequency)
    nRad = study.Ns
    nCoeff = len(study.p_mesh[0, 0].coefficients)
    nCoeffV = len(study.u_mesh[0, 0].coefficients)
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    velocityMesh = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    velocityMesh_Y = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    for freq in range(Nfft):
        for rad in range(nRad):
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
            velocityMesh[freq, rad, :] = study.u_mesh[freq, rad].coefficients
            velocityMesh_Y[freq, rad, :] = study.u_mesh_Y[freq, rad].coefficients
    return pressureMesh, velocityMesh, velocityMesh_Y


def storePressureMeshResults(acoustic_study):
    study = acoustic_study
    Nfft = len(study.frequency)
    nRad = study.Ns
    nCoeff = len(study.p_mesh[0, 0].coefficients)
    nCoeffV = len(study.u_mesh[0, 0].coefficients)
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    velocityMesh = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    for freq in range(Nfft):
        for rad in range(nRad):
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
            velocityMesh[freq, rad, :] = study.u_mesh[freq, rad].coefficients
    return pressureMesh, velocityMesh


# LOAD PRESSURE
def loadPressureMeshResults(obj, pressureMesh):
    Nfft  = np.shape(pressureMesh)[0]
    nRad  = np.shape(pressureMesh)[1]

    for f in range(Nfft):
        for rs in range(nRad):
            obj.p_mesh[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, coefficients=pressureMesh[f, rs, :])
            dofCount = obj.spaceU_freq[rs].grid_dof_count
            coeff_radSurf = np.ones(dofCount, dtype=complex) * obj.coeff_radSurf[f, rs, :dofCount]
            spaceU = bempp_cl.api.function_space(obj.grid_sim, "DP", 0,
                                              segments=[obj.radiatingElement[rs]])
            u_total = bempp_cl.api.GridFunction(spaceU, coefficients=-coeff_radSurf)
            obj.u_mesh[f, rs] = u_total

        obj.p_total_mesh[f] = bempp_cl.api.GridFunction(obj.spaceP,
                                                     coefficients=np.sum(pressureMesh[f, :, :], 0))
    return None


def loadPointSourceBEM(obj, pressureMesh, velocityMesh, velocityMesh_Y):
    Nfft  = np.shape(pressureMesh)[0]
    nRad  = np.shape(pressureMesh)[1]
    for f in range(Nfft):
        for rs in range(nRad):
            obj.p_mesh[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                       coefficients=pressureMesh[f, rs, :])
            obj.u_mesh[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                       coefficients=velocityMesh[f, rs, :])        
            if obj.admittanceCoeff is not None:
                obj.u_mesh_Y[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                             coefficients=velocityMesh_Y[f, rs, :])
    return None
