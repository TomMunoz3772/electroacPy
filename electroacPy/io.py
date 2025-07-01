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
             frequency   = frequency,
             driver      = sim.driver,
             laser       = sim.vibrometry,
             crossover   = sim.crossover,
             enclosure   = sim.enclosure,
             radiator_id = sim.radiator_id,
             c           = sim.c,
             rho         = sim.rho)

    # Save evaluations / bem
    for study in sim.acoustic_study:
        if bool(sim.evaluation) is True:
            obsTmp = sim.evaluation[study]
            np.savez(join(projectPath, 'evs_{}'.format(study)),
                     frequency            = obsTmp.frequency,
                     setup                = obsTmp.setup,
                     referenceStudy       = obsTmp.referenceStudy,)
        else:
            pass

        # BEM studies
        studyTmp = sim.acoustic_study[study]
        
        # prepare to store mesh pressure and velocity 
        # store PS with admittance values
        if studyTmp.admittanceCoeff is not None and hasattr(studyTmp, "xSource"):
            pressureArrayAcs, velocityArrayAcs, velocityArrayAcs_Y, propagCoeff, Yn_c = storePressureMeshResults_PSADM(studyTmp)
        # store PS without admittance values
        elif studyTmp.admittanceCoeff is None and hasattr(studyTmp, "xSource"):
            pressureArrayAcs, velocityArrayAcs, propagCoeff = storePressureMeshResults_PS(studyTmp)
            velocityArrayAcs_Y = None
            Yn_c = None
        # store BEM with admittance values
        elif studyTmp.admittanceCoeff is not None:
            pressureArrayAcs, velocityArrayAcs, Yn_c = storePressureMeshResults_ADM(studyTmp)
        # store BEM without admittance values
        else:
            pressureArrayAcs, velocityArrayAcs = storePressureMeshResults(studyTmp)
            Yn_c = None            
            
        # do a copy of mesh in save folder
        if studyTmp.meshPath is not None:
            copy2(studyTmp.meshPath, projectPath)
            mesh_filename = os.path.basename(studyTmp.meshPath)
        else:
            mesh_filename = None
            
        # EXTERIOR
        if hasattr(studyTmp, "xSource"):
            np.savez(join(projectPath, 'acs_{}'.format(study)),
                     meshPath           = studyTmp.meshPath,
                     mesh_filename      = mesh_filename,
                     radiatingElement   = studyTmp.radiatingElement,
                     impedanceElement   = studyTmp.impedanceSurfaceIndex ,
                     QSource            = studyTmp.QSource,
                     xSource            = studyTmp.xSource,
                     isComputed         = studyTmp.isComputed,
                     vertices           = studyTmp.vertices,
                     domain             = studyTmp.domain,
                     kwargs             = studyTmp.kwargs,
                     pressureArrayAcs   = pressureArrayAcs,
                     velocityArrayAcs   = velocityArrayAcs,
                     velocityArrayAcs_Y = velocityArrayAcs_Y,
                     propagCoeff        = propagCoeff,
                     Yn_c               = Yn_c, 
                     LEM_enclosures     = studyTmp.LEM_enclosures,
                     radiator           = studyTmp.radiator,
                     c_0                = studyTmp.c_0,
                     rho_0              = studyTmp.rho_0
                     )
        else:
            np.savez(join(projectPath, 'acs_{}'.format(study)),
                     meshPath             = studyTmp.meshPath,
                     mesh_filename        = mesh_filename,
                     radiatingElement     = studyTmp.radiatingElement,
                     velocity             = array(studyTmp.velocity, dtype=object),
                     isComputed           = studyTmp.isComputed,
                     coeff_radSurf        = studyTmp.coeff_radSurf,
                     vertices             = studyTmp.vertices,
                     domain               = studyTmp.domain,
                     kwargs               = studyTmp.kwargs,
                     pressureArrayAcs     = pressureArrayAcs,
                     Yn_c                 = Yn_c,
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
    LS: loudspeakerSystem object,
        Study, evaluations and LEM setup.
    """

    dataLEM = np.load(join(pathToProject, 'LEM.npz'), allow_pickle=True)

    # create loudspeaker system
    LS             = loudspeakerSystem(dataLEM['frequency'])
    LS.driver      = dataLEM['driver'].item()
    LS.vibrometry  = dataLEM['laser'].item()
    LS.enclosure   = dataLEM['enclosure'].item()
    LS.crossover   = dataLEM['crossover'].item()
    LS.radiator_id = dataLEM['radiator_id'].item()
    
    try: # to keep it compatible with previous datasets
        LS.c   = dataLEM['c']
        LS.rho = dataLEM['rho']
    except:
        LS.c   = 343
        LS.rho = 1.22

    # import studies and evaluations
    file_list  = os.listdir(pathToProject)
    acs_files  = [file for file in file_list if file.startswith('acs')]
    study_name = []
    for i in range(len(acs_files)):
        study_name.append(acs_files[i][4:-4])

    for study in study_name:
        # load acoustic_study
        data_acs     = np.load(join(pathToProject, 'acs_{}.npz'.format(study)),
                               allow_pickle=True)
        meshName     = data_acs['mesh_filename'].item()
        isComputed   = data_acs['isComputed'].item()
        domain       = data_acs["domain"].item()
        kwargs       = data_acs['kwargs'].item()

        try:
            enclosures   = list(data_acs['LEM_enclosures'])
        except:
            enclosures    = {}
        try:
            radiator   = list(data_acs['radiator'])
        except:
            radiator    = {}
    
        if "xSource" in data_acs and meshName is not None: # pointSourceBEM
            xSource = data_acs["xSource"]
            QSource = data_acs["QSource"]
            physics_acs = pointSourceBEM(join(pathToProject, meshName),
                                         xSource, QSource, 
                                         LS.frequency, domain=domain, 
                                         c_0=LS.c, rho_0=LS.rho, 
                                         **kwargs)
            physics_acs.isComputed     = isComputed
            physics_acs.LEM_enclosures = enclosures
            physics_acs.radiator       = radiator
            loadPointSourceBEM(physics_acs, data_acs["pressureArrayAcs"],
                               data_acs["velocityArrayAcs"], 
                               data_acs["velocityArrayAcs_Y"], 
                               data_acs["propagCoeff"], 
                               data_acs["Yn_c"])
            LS.acoustic_study[study] = physics_acs
            
        elif "xSource" in data_acs and meshName is None: # pointSource
            xSource = data_acs["xSource"]
            QSource = data_acs["QSource"]
            physics_acs = pointSource(xSource, QSource, LS.frequency,
                                      c=LS.c, rho=LS.rho,
                                      boundary_conditions=kwargs["boundary_conditions"],
                                      **kwargs)
            physics_acs.isComputed = isComputed
            physics_acs.LEM_enclosures = enclosures
            physics_acs.radiator       = radiator
            LS.acoustic_study[study] = physics_acs
        else: # BEM
            radSurf      = data_acs['radiatingElement']
            surfVelocity = data_acs['velocity']  
            physics_acs = bem(join(pathToProject, meshName),
                                  radSurf,
                                  surfVelocity,
                                  LS.frequency,
                                  domain=domain,
                                  c_0=LS.c, rho_0=LS.rho,
                                  **kwargs)
            physics_acs.isComputed     = isComputed
            physics_acs.LEM_enclosures = enclosures
            physics_acs.radiator       = radiator
            loadPressureMeshResults(physics_acs, data_acs['pressureArrayAcs'],
                                    data_acs["Yn_c"])
            LS.acoustic_study[study] = physics_acs

        # load evaluations
        try:
            try:
                data_evs = np.load(join(pathToProject, 'evs_{}.npz'.format(study)), allow_pickle=True)
            except:
                data_evs = np.load(join(pathToProject, 'obs_{}.npz'.format(study)), allow_pickle=True)
                print("Legacy evs loaded")
            physics_evs = evs(physics_acs)
            physics_evs.setup = data_evs["setup"].item()
            physics_evs.frequency = data_evs["frequency"]
            physics_evs.referenceStudy = data_evs["referenceStudy"].item()
            LS.evaluation[study] = physics_evs
        except:
            print('No evaluation to load')
    return LS


#%% HELPERS
def storePressureMeshResults_PSADM(acoustic_study):
    study       = acoustic_study
    Nfft        = len(study.frequency)
    nRad        = study.Ns
    nCoeff      = len(study.p_mesh[0, 0].coefficients)
    nCoeffV     = len(study.u_mesh[0, 0].coefficients)
    propagCoeff = np.zeros((Nfft, nRad, study.spaceP.grid_dof_count), dtype=complex)
    Yn_c        = np.zeros((Nfft, study.spaceP.grid_dof_count), dtype=complex)

    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    velocityMesh = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    velocityMesh_Y = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    for freq in range(Nfft):
        Yn_c[freq, :] = study.Yn_c[freq]
        for rad in range(nRad):
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
            velocityMesh[freq, rad, :] = study.u_mesh[freq, rad].coefficients
            velocityMesh_Y[freq, rad, :] = study.u_mesh_Y[freq, rad].coefficients
            propagCoeff[freq, rad, :] = study.propag_function[freq, rad].coefficients
    return pressureMesh, velocityMesh, velocityMesh_Y, propagCoeff, Yn_c


def storePressureMeshResults_PS(acoustic_study):
    study = acoustic_study
    Nfft = len(study.frequency)
    nRad = study.Ns
    nCoeff = len(study.p_mesh[0, 0].coefficients)
    nCoeffV = len(study.u_mesh[0, 0].coefficients)
    propagCoeff = np.zeros((Nfft, nRad, nCoeff), dtype=complex)
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    velocityMesh = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    velocityMesh_Y = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    for freq in range(Nfft):
        for rad in range(nRad):
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
            velocityMesh[freq, rad, :] = study.u_mesh[freq, rad].coefficients
            # velocityMesh_Y[freq, rad, :] = study.u_mesh_Y[freq, rad].coefficients
            propagCoeff[freq, rad, :] = study.propag_function[freq, rad].coefficients
    return pressureMesh, velocityMesh, propagCoeff

def storePressureMeshResults(acoustic_study):
    study = acoustic_study
    Nfft = len(study.frequency)
    nRad = study.Ns
    nCoeff = len(study.p_mesh[0, 0].coefficients)
    
    # get max u_mesh len (because acoustic radiators can be of different size, i.e. mesh elements)
    radSurf_elements = []
    for i in range(np.shape(study.u_mesh)[1]):
        tmp_length = len(study.u_mesh[0, i].coefficients)
        radSurf_elements.append(tmp_length)
    nCoeffV = np.max(radSurf_elements)
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    velocityMesh = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    for freq in range(Nfft):
        for rad in range(nRad):
            tmp_length = radSurf_elements[rad]
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
            velocityMesh[freq, rad, :tmp_length] = study.u_mesh[freq, rad].coefficients
    return pressureMesh, velocityMesh

def storePressureMeshResults_ADM(acoustic_study):
    study = acoustic_study
    Nfft = len(study.frequency)
    nRad = study.Ns
    nCoeff = len(study.p_mesh[0, 0].coefficients)
    # nCoeffV = len(study.u_mesh[0, 0].coefficients)
    Yn_c   = np.zeros((Nfft, nCoeff), dtype=complex)
    
    # get max u_mesh len (because acoustic radiators can be of different size, i.e. mesh elements)
    radSurf_elements = []
    for i in range(np.shape(study.u_mesh)[1]):
        tmp_length = len(study.u_mesh[0, i].coefficients)
        radSurf_elements.append(tmp_length)
    nCoeffV = np.max(radSurf_elements)
    
    # store pressure
    pressureMesh = np.zeros([Nfft, nRad, nCoeff], dtype=complex)
    velocityMesh = np.zeros([Nfft, nRad, nCoeffV], dtype=complex)
    for freq in range(Nfft):
        Yn_c[freq, :] = study.Yn_c[freq]
        for rad in range(nRad):
            tmp_length = radSurf_elements[rad]
            pressureMesh[freq, rad, :] = study.p_mesh[freq, rad].coefficients
            velocityMesh[freq, rad, :tmp_length] = study.u_mesh[freq, rad].coefficients
    return pressureMesh, velocityMesh, Yn_c

# LOAD PRESSURE
def loadPressureMeshResults(obj, pressureMesh, Yn_c):
    Nfft  = np.shape(pressureMesh)[0]
    nRad  = np.shape(pressureMesh)[1]
    obj.Yn_c = Yn_c
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


def loadPointSourceBEM(obj, pressureMesh, velocityMesh, velocityMesh_Y, propagCoeff, Yn_c):
    Nfft  = np.shape(pressureMesh)[0]
    nRad  = np.shape(pressureMesh)[1]
    obj.Yn_c = Yn_c
    for f in range(Nfft):
        for rs in range(nRad):
            obj.p_mesh[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                       coefficients=pressureMesh[f, rs, :])
            obj.u_mesh[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                       coefficients=velocityMesh[f, rs, :])
            obj.propag_function[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                                   coefficients=propagCoeff[f, rs, :])
            if obj.admittanceCoeff is not None:
                obj.u_mesh_Y[f, rs] = bempp_cl.api.GridFunction(obj.spaceP, 
                                                             coefficients=velocityMesh_Y[f, rs, :])
    return None
