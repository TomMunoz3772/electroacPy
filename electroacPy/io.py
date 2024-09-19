"""
Save and Load projects
"""
import numpy as np
import os
from shutil import copy2
from os.path import join
from numpy import asanyarray as array
from electroacPy import loudspeakerSystem
from electroacPy.acousticSim.exteriorAcoustics.bem import bem as bem_ext
from electroacPy.acousticSim.exteriorAcoustics.observations import observations as  obs_ext
from electroacPy.acousticSim.interiorAcoustics.bem import bem as bem_int
from electroacPy.acousticSim.interiorAcoustics.observations import observations as  obs_int

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
        try:
            obsTmp = sim.observation[study]
            pressureArrayObs, nMic = storePressureMicResults(obsTmp)
            np.savez(join(projectPath, 'obs_{}'.format(study)),
                     observationName      = obsTmp.observationName,
                     observationType      = obsTmp.observationType,
                     computedObservations = obsTmp.computedObservations,
                     xMic                 = array(obsTmp.xMic, dtype=object),
                     labels               = array(obsTmp.labels, dtype=object),
                     L                    = array(obsTmp.L, dtype=object),
                     W                    = array(obsTmp.W, dtype=object),
                     Lx                   = array(obsTmp.Lx, dtype=object),
                     Ly                   = array(obsTmp.Ly, dtype=object),
                     Lz                   = array(obsTmp.Lz, dtype=object),
                     planarPlane          = obsTmp.planarPlane,
                     add2Plotter          = obsTmp.add2Plotter,
                     theta                = array(obsTmp.theta, dtype=object),
                     polarPlane           = obsTmp.polarPlane,
                     referenceStudy       = obsTmp.referenceStudy,
                     pressureArrayObs     = pressureArrayObs,
                     nMic                 = nMic
                     )
        except:
            None

        # BEM studies
        studyTmp = sim.acoustic_study[study]
        pressureArrayAcs = storePressureMeshResults(studyTmp)
        copy2(studyTmp.meshPath, projectPath)
        mesh_filename = os.path.basename(studyTmp.meshPath)
        # EXTERIOR
        np.savez(join(projectPath, 'acs_{}'.format(study)),
                 meshPath            = studyTmp.meshPath,
                 mesh_filename       = mesh_filename,
                 radiatingSurface    = studyTmp.radiatingSurface,
                 surfaceVelocity     = array(studyTmp.surfaceVelocity, dtype=object),
                 radiation_direction = array(studyTmp.radiation_direction, dtype=object),
                 vibrometry_points   = array(studyTmp.vibrometry_points, dtype=object),
                 nRad_S              = studyTmp.nRad_S,
                 isComputed          = studyTmp.isComputed,
                 coeff_radSurf       = studyTmp.coeff_radSurf,
                 vertices            = studyTmp.vertices,
                 boundary            = studyTmp.boundary,
                 offset              = studyTmp.offset,
                 pressureArrayAcs    = pressureArrayAcs,
                 identifier          = studyTmp.identifier,
                 LEM_enclosures      = studyTmp.LEM_enclosures,
                 radiator            = studyTmp.radiator,
                 c_0                 = studyTmp.c_0,
                 rho_0               = studyTmp.rho_0
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
        data_acs   = np.load(join(pathToProject, 'acs_{}.npz'.format(study)), allow_pickle=True)
        study_ID = data_acs['identifier'].item()
        meshName   = data_acs['mesh_filename'].item()
        if study_ID == 'BEM_EXT':
            radSurf      = data_acs['radiatingSurface']
            surfVelocity = data_acs['surfaceVelocity']
            try:
                boundary = data_acs['boundary'].item()
            except:
                boundary = list(data_acs['boundary'])
            try:
                offset = list(data_acs['offset'])  # might cause issues ['n', 'o', 'n', 'e']
                # print(offset)
            except:
                offset = data_acs['offset'].item()
                # print(offset)
            isComputed   = data_acs['isComputed'].item()
            try:
                enclosures   = list(data_acs['LEM_enclosures'])
            except:
                enclosures    = {}
            
            try:
                radiator   = list(data_acs['radiator'])
            except:
                radiator    = {}

            # vibrometry_points = list(np.squeeze(data_acs['vibrometry_points']))
            try:
                # vibrometry_points = np.squeeze(data_acs['vibrometry_points'])
                vibrometry_points = data_acs['vibrometry_points']
            except:
                vibrometry_points = False

            # print(vibrometry_points[0].shape)

            try:
                radiation_direction = data_acs['radiation_direction'].item()
            except:
                radiation_direction = list(data_acs['radiation_direction'])
            
            physics_acs = bem_ext(join(pathToProject, meshName),
                                  radSurf,
                                  surfVelocity,
                                  LS.frequency,
                                  vibrometry_points,
                                  radiation_direction,
                                  boundary,
                                  offset, c_0=LS.c, rho_0=LS.rho)
            physics_acs.isComputed     = isComputed
            physics_acs.LEM_enclosures = enclosures
            physics_acs.radiator       = radiator
            loadPressureMeshResults(physics_acs, data_acs['pressureArrayAcs'])
            LS.acoustic_study[study] = physics_acs

            # load observations
            try:
                data_obs = np.load(join(pathToProject, 'obs_{}.npz'.format(study)), allow_pickle=True)
                physics_obs = obs_ext(physics_acs)
                variable_names = ['observationName', 'observationType', 'computedObservations',
                  'xMic', 'labels', 'L', 'W', 'Lx', 'Ly', 'Lz', 'planarPlane', 'add2Plotter',
                  'theta', 'polarPlane', 'referenceStudy', 'pressureArrayObs',
                  'nMic']
                for var_name in variable_names:
                    try:
                        saved_value = list(data_obs[var_name])
                        setattr(physics_obs, var_name, saved_value)
                    except:
                        try:
                            setattr(physics_obs, var_name, ['N/A'] * len(list(data_obs['observationName'])))
                            print(f"{var_name}: initialized variable")
                        except:
                            print('No observation to load - check if observationName is empty')
                loadPressureMicResults(physics_obs, data_obs['pressureArrayObs'], data_obs['nMic'])
                LS.observation[study] = physics_obs
            except:
                print('No observation to load')


        if study_ID == 'BEM_INT':
            radSurf = data_acs['radiatingSurface']
            surfVelocity = data_acs['surfaceVelocity']
            absorbingSurface = data_acs['absorbingSurface']
            surfaceImpedance = data_acs['surfaceImpedance']
            isComputed = data_acs['isComputed'].item()
            try:
                enclosures = list(data_acs['LEM_enclosures'])
            except:
                enclosures = {}

            try:
                radiator = list(data_acs['radiator'])
            except:
                radiator = {}

            try:
                vibrometry_points = list(np.squeeze(data_acs['vibrometry_points']))
            except:
                vibrometry_points = False

            try:
                radiation_direction = data_acs['radiation_direction'].item()
            except:
                radiation_direction = list(data_acs['radiation_direction'])

            physics_acs = bem_int(join(pathToProject, meshName),
                                  radSurf,
                                  surfVelocity,
                                  LS.frequency,
                                  absorbingSurface,
                                  surfaceImpedance,
                                  vibrometry_points,
                                  radiation_direction, c_0=LS.c, rho_0=LS.rho)
            physics_acs.isComputed = isComputed
            physics_acs.LEM_enclosures = enclosures
            physics_acs.radiator = radiator
            loadPressureMeshResults(physics_acs, data_acs['pressureArrayAcs'])
            LS.acoustic_study[study] = physics_acs

            # load observations
            try:
                data_obs = np.load(join(pathToProject, 'obs_{}.npz'.format(study)), allow_pickle=True)
                physics_obs = obs_int(physics_acs)   ## OBJECT REFERS AS OBS_BEM_EXT ????
                physics_obs.observationName = list(data_obs['observationName'])
                physics_obs.observationType = list(data_obs['observationType'])
                physics_obs.computedObservations = list(data_obs['computedObservations'])
                physics_obs.xMic = list(data_obs['xMic'])
                physics_obs.labels = list(data_obs['labels'])
                physics_obs.L = list(data_obs['L'])
                physics_obs.W = list(data_obs['W'])
                physics_obs.planarPlane = list(data_obs['planarPlane'])
                physics_obs.add2Plotter = list(data_obs['add2Plotter'])
                physics_obs.theta = list(data_obs['theta'])
                physics_obs.polarPlane = list(data_obs['polarPlane'])
                physics_obs.referenceStudy = data_obs['referenceStudy'].item()
                loadPressureMicResults(physics_obs, data_obs['pressureArrayObs'], data_obs['nMic'])
                LS.observation[study] = physics_obs
            except:
                print('No observation to load')

    return LS


# def storeMicPosition(xMic):
#     nObs = len(xMic)
#     micLen = []
#     for i in range(nObs):
#         micLen.append(np.max(xMic[i].shape))
#     maxMic = np.max(micLen)
#     micToStore = np.zeros([nObs, maxMic, 3])
#     for i in range(nObs):
#         micToStore[i, :micLen[i], :] = xMic[i][:micLen[i], :]
#     print(np.shape(micToStore))
#     return micToStore
#
# def loadMicPosition(storedMic, nMic):
#     xMic = []
#     for i in range(storedMic.shape[0]):
#         xMic.append(storedMic)
#     return None


def storePressureMicResults(observation):
    obs = observation
    nObs = len(obs.observationName)       # number of observtations
    Nfft = len(obs.bemObject.freq_array)  # number of frequency bins
    nRad = obs.bemObject.nRad_S      # number of radiating surfaces
    nMic = np.zeros(nObs)                 # number of microphone for each radSurf
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
    Nfft = len(study.freq_array)
    nRad = study.nRad_S
    nVert = study.vertices
    nCoeff = study.spaceP.grid_dof_count
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    for freq in range(Nfft):
        for rad in range(nRad):
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
    
    # store velocity attribution
    # radSurf_nVert = np.zeros(nRad, dtype=complex)  # will store number of triangle part of each radiating surface
    # for rad in range(nRad):
    #     radSurf_nVert[rad] = len(study.u_total_list_array[0, rad].coefficients)
    # maxVert = np.max(radSurf_nVert) # what is the maximum number of triangles?
    # velocityMesh = np.zeros([Nfft, nRad, maxVert])
    # for freq in range(Nfft):
    #     for rad in range(nRad):
    #         velocityMesh[freq, rad, :radSurf_nVert[rad]] = study.
    return pressureMesh

def loadPressureMeshResults(obj, pressureMesh):
    Nfft  = np.shape(pressureMesh)[0]
    nRad  = np.shape(pressureMesh)[1]
    nVert = np.shape(pressureMesh)[2]

    for f in range(Nfft):
        for rs in range(nRad):
            obj.p_mesh[f, rs] = bempp.api.GridFunction(obj.spaceP, coefficients=pressureMesh[f, rs, :])
            dofCount = obj.spaceU_list[rs].grid_dof_count
            coeff_radSurf = np.ones(dofCount, dtype=complex) * obj.coeff_radSurf[f, rs, :dofCount]
            spaceU = bempp.api.function_space(obj.grid, "DP", 0,
                                              segments=[obj.radiatingSurface[rs]])
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
