"""
Save and Load projects
"""
import numpy as np
import os
from shutil import copy2
from os.path import join
from numpy import asanyarray as array
from electroacPy import LS2 as loudspeakerSystem
from electroacPy.acousticSim.bem import bem
from electroacPy.acousticSim.observations import observations as  obs
import bempp.api

def save(projectPath, loudspeakerSystem):
    """
    Save loudspeaker system simulation in project folder.

    :param loudspeakerSystem:
    :return:
    *.npz archive
    """

    if not os.path.exists(projectPath):
        os.mkdir(projectPath)

    sim = loudspeakerSystem

    # global variables
    frequency = sim.frequency

    # Save lumped element / analytic simulations
    np.savez(join(projectPath, 'LEM'),
             frequency=frequency,
             driver=sim.driver,
             laser = sim.laser_acc,
             crossover=sim.crossover,
             enclosure=sim.enclosure,
             radiator_id=sim.radiator_id,
             c=sim.c,
             rho=sim.rho)

    # Save observations / bem
    for study in sim.acoustic_study:
        if bool(sim.observation) is True:
            obsTmp = sim.observation[study]
            # pressureArrayObs, nMic = storePressureMicResults(obsTmp)
            np.savez(join(projectPath, 'obs_{}'.format(study)),
                     frequency            = obsTmp.frequency,
                     setup                = obsTmp.setup,
                     referenceStudy       = obsTmp.referenceStudy,)
        else:
            pass

        # BEM studies
        studyTmp = sim.acoustic_study[study]
        pressureArrayAcs = storePressureMeshResults(studyTmp)
        copy2(studyTmp.meshPath, projectPath)
        mesh_filename = os.path.basename(studyTmp.meshPath)
        # EXTERIOR
        np.savez(join(projectPath, 'acs_{}'.format(study)),
                 meshPath             = studyTmp.meshPath,
                 mesh_filename        = mesh_filename,
                 radiatingElement     = studyTmp.radiatingElement,
                 velocity             = array(studyTmp.velocity, dtype=object),
                 isComputed           = studyTmp.isComputed,
                 coeff_radSurf        = studyTmp.coeff_radSurf,
                 vertices             = studyTmp.vertices,
                 kwargs               = studyTmp.kwargs,
                 pressureArrayAcs     = pressureArrayAcs,
                 LEM_enclosures       = studyTmp.LEM_enclosures,
                 radiator             = studyTmp.radiator,
                 c_0                  = studyTmp.c_0,
                 rho_0                = studyTmp.rho_0
                 )
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
    LS: loudspeakerSystem object
    """

    dataLEM = np.load(join(pathToProject, 'LEM.npz'), allow_pickle=True)

    # create loudspeaker system
    LS             = loudspeakerSystem(dataLEM['frequency'])
    LS.driver      = dataLEM['driver'].item()
    LS.laser_acc   = dataLEM['laser'].item()
    LS.enclosure   = dataLEM['enclosure'].item()
    LS.crossover   = dataLEM['crossover'].item()
    LS.radiator_id = dataLEM['radiator_id'].item()
    try: # to keep it compatible with previous datasets
        LS.c   = dataLEM['c']
        LS.rho = dataLEM['rho']
    except:
        LS.c   = 343
        LS.rho = 1.22

    # import studies and observations
    file_list  = os.listdir(pathToProject)
    acs_files  = [file for file in file_list if file.startswith('acs')]
    obs_files  = [file for file in file_list if file.startswith('obs')]
    study_name = []
    for i in range(len(acs_files)):
        study_name.append(acs_files[i][4:-4])

    for study in study_name:
        # load acoustic_study
        data_acs   = np.load(join(pathToProject, 'acs_{}.npz'.format(study)),
                             allow_pickle=True)
        meshName   = data_acs['mesh_filename'].item()
        radSurf      = data_acs['radiatingElement']
        surfVelocity = data_acs['velocity']
        isComputed   = data_acs['isComputed'].item()
        kwargs       = data_acs['kwargs'].item()
        try:
            enclosures   = list(data_acs['LEM_enclosures'])
        except:
            enclosures    = {}
        try:
            radiator   = list(data_acs['radiator'])
        except:
            radiator    = {}
        physics_acs = bem(join(pathToProject, meshName),
                              radSurf,
                              surfVelocity,
                              LS.frequency,
                              **kwargs)
        physics_acs.isComputed     = isComputed
        physics_acs.LEM_enclosures = enclosures
        physics_acs.radiator       = radiator
        loadPressureMeshResults(physics_acs, data_acs['pressureArrayAcs'])
        LS.acoustic_study[study] = physics_acs

        # load observations
        try:
            data_obs = np.load(join(pathToProject, 'obs_{}.npz'.format(study)), allow_pickle=True)
            physics_obs = obs(physics_acs)
            physics_obs.setup = data_obs["setup"].item()
            physics_obs.frequency = data_obs["frequency"]
            physics_obs.referenceStudy = data_obs["referenceStudy"].item()
            LS.observation[study] = physics_obs
        except:
            print('No observation to load')
    return LS


def storePressureMicResults(observation):
    obs = observation
    nObs = len(obs.setup)
    Nfft = len(obs.bemObject.frequency)  # number of frequency bins
    nRad = obs.bemObject.Ns              # number of radiating surfaces
    nMic = np.zeros(nObs)                # number of microphone for each radSurf
    for i in range(nObs):
        nMic[i] = obs.xMic[i].shape[0]
    maxMic = int(np.max(nMic))

    pressureMicrophone = np.zeros([nObs, Nfft, maxMic, nRad], dtype=complex)
    for obsTmp in range(nObs):
        for freq in range(Nfft):
            for rad in range(nRad):
                pressureMicrophone[obsTmp, freq, :int(nMic[obsTmp]), rad] = obs.pMicArray[obsTmp][freq, :, rad]
    return pressureMicrophone, nMic

def storePressureMeshResults(acoustic_study):
    study = acoustic_study
    Nfft = len(study.frequency)
    nRad = study.Ns
    nVert = study.vertices
    nCoeff = study.spaceP.grid_dof_count
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    for freq in range(Nfft):
        for rad in range(nRad):
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
    return pressureMesh

def loadPressureMeshResults(obj, pressureMesh):
    Nfft  = np.shape(pressureMesh)[0]
    nRad  = np.shape(pressureMesh)[1]
    nVert = np.shape(pressureMesh)[2]

    for f in range(Nfft):
        for rs in range(nRad):
            obj.p_mesh[f, rs] = bempp.api.GridFunction(obj.spaceP, coefficients=pressureMesh[f, rs, :])
            dofCount = obj.spaceU_freq[rs].grid_dof_count
            coeff_radSurf = np.ones(dofCount, dtype=complex) * obj.coeff_radSurf[f, rs, :dofCount]
            spaceU = bempp.api.function_space(obj.grid_sim, "DP", 0,
                                              segments=[obj.radiatingElement[rs]])
            u_total = bempp.api.GridFunction(spaceU, coefficients=-coeff_radSurf)
            obj.u_mesh[f, rs] = u_total

        # print('SHAPE: ', np.sum(pressureMesh[f, :, :], 0).shape)
        obj.p_total_mesh[f] = bempp.api.GridFunction(obj.spaceP,
                                                     coefficients=np.sum(pressureMesh[f, :, :], 0))

    return None


def loadPressureMicResults(obj, pressureMic, nMic):
    nObs   = np.shape(pressureMic)[0]
    Nfft   = np.shape(pressureMic)[1]
    maxMic = np.shape(pressureMic)[2]
    nRad   = np.shape(pressureMic)[3]
    pMic = []
    pMicArray = []
    for i in range(nObs):
        pMicArray.append(pressureMic[i, :, :int(nMic[i]), :])
        pMic.append(np.sum(pressureMic[i, :, :int(nMic[i]), :], 2))
    obj.pMic = pMic
    obj.pMicArray = pMicArray
    return None
