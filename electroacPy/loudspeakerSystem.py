#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 16:03:02 2023

@author: tom.munoz
"""

#==================================================
#
#    Loudspeaker system simulation class
#
#==================================================
# BEM acoustic sim
from electroacPy.acousticSim.bem import bem
from electroacPy.acousticSim.evaluations import evaluations as evs_bem
from electroacPy.acousticSim.evaluations import getPressure
from electroacPy.acousticSim.postProcessing import postProcess as pp

# PS acoustic sim
from electroacPy.acousticSim.pointSource import pointSource, pointSourceBEM

# Lumped element
from electroacPy.speakerSim.electroAcousticDriver import electroAcousticDriver, loadLPM
from electroacPy.speakerSim.enclosureDesign import speakerBox
from electroacPy.speakerSim.filterDesign import xovers

# Vibrometry
from electroacPy.measurements.laserVibrometry import laserVibrometry_UFF as laser_v_uff
from electroacPy.measurements.laserVibrometry import laserVibrometry as laser_v

# general
from electroacPy.global_ import air
from electroacPy.general.freqop import freq_log10
from electroacPy.general.gain import dB
from electroacPy.general.plot import bempp_grid_mesh

# external libraries
import numpy as np
import matplotlib.pyplot as plt


## CLASS
class loudspeakerSystem:
    """
    Automate the use of acousticSim + speakerSim modules specifically for loudspeaker
    design.
    
    Parameters
    ----------
    frequencyRange: numpy array, optional
        Delimitation of the study. Default is freq_log10(20, 2500, 50).
    
    c: float, optional
        Speed of sound in the propagation medium. Default is 343 m/s.
    
    rho: float, optional
        Medium density. Default if 1.22 kg/m^3
        
    Returns
    -------
    None
    """ 
    def __init__(self, frequencyRange=freq_log10(20, 2500, 50),
                 c=air.c, rho=air.rho, **kwargs):
        # global variables
        self.frequency = frequencyRange
        self.c = c
        self.rho = rho

        # simulations
        self.driver         = {}
        self.enclosure      = {}
        self.vibrometry     = {}
        self.crossover      = {}
        self.acoustic_study = {}
        self.evaluation     = {}
        self.results        = {}

        # help
        self.radiator_id    = {}

        if "compiler_output" in kwargs:
            if kwargs["compiler_output"] is True:
                import os
                os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
            else:
                pass

    ## ===============================
    # %% LUMPED ELEMENTS / LOUDSPEAKER
    def lem_driver(self, name, U, Le, Re, Cms, Mms, Rms, Bl, Sd, ref2bem=None):
        """
        Add electro-dynamic driver from its Thiele/Small parameters. Can act as a radiator.
        :param name:
        :param U:
        :param Le:
        :param Re:
        :param Cms:
        :param Mms:
        :param Rms:
        :param Bl:
        :param Sd:
        :param ref2bem:
        :return:
        """
        physics = electroAcousticDriver(U, Le, Re, Cms, Mms, Rms, Bl, Sd,
                                        self.frequency, self.c, self.rho)
        physics.ref2bem = ref2bem
        self.driver[name] = physics
        self.radiator_id[name] = 'EAD'
        return None

    def lem_driverImport(self, name, lpm_data, U=1, ref2bem=None, 
                           LeZero=False):
        """
        Add electro-dynamic driver from a LPM file (Thiele/Small parameters from Klippel bench). Can act as a radiator
        :param name:
        :param lpm_data:
        :param U:
        :param ref2bem:
        :return:
        """
        physics = loadLPM(lpm_data, self.frequency, U=U, c=self.c, rho=self.rho, LeZero=LeZero)
        physics.ref2bem = ref2bem
        self.driver[name] = physics
        self.radiator_id[name] = 'EAD'
        return None

    def lem_velocity(self, name, ref2bem, U=1, v=1):
        physics = electroAcousticDriver(U, 1, 1, 1, 1, 1, 1, 1,
                                        self.frequency, self.c, self.rho) # "false" loudspeaker
        physics.ref2bem = ref2bem
        physics.v = np.ones(len(self.frequency), dtype=complex) * v
        physics.Q = np.ones(len(self.frequency), dtype=complex) * v
        self.driver[name] = physics
        self.radiator_id[name] = 'EAD'
        return None

    def lem_enclosure(self, name, Vb, Qab=120, Qal=30, setDriver=None, Nd=1, 
                      wiring="parallel", ref2bem=None, **kwargs):
        
        physics = speakerBox(Vb, frequencyRange=self.frequency, 
                             c=self.c, rho=self.rho,
                             Qab=Qab, Qal=Qal, **kwargs)
        physics.ref2bem = ref2bem
        self.enclosure[name] = physics
        self.radiator_id[name] = 'SPKBOX'

        if setDriver != None:
            physics.getDriverResponse(self.driver[setDriver], 
                                      Nd=Nd, wiring=wiring)
            physics.whichDriver = setDriver
        return None

    def vibrometry_data(self, name, file_path, rotation=[0, 0, 0],
                        ref2bem=None, useAverage=False, inputVoltage=1):
        """
        Add acceleration data to the study. Meant to be used as a radiator - ref2bem strongly recommended.

        Parameters
        ----------
        name : str
            Reference.
        file_path : str
            Path to *.uff vibrometry data.
        rotation : list of float, optional
            Rotate measured data to match mesh direction. The default is [0, 0, 0].
        ref2bem : int, optional
            Reference to mesh physical group. The default is None.
        useAverage : Bool, optional
            Uses the average instead of individual element acceleration. 
            The default is False.
        inputVoltage : float, optional
            Optional scaling. By default, electroacPy uses the transfer function of the 
            acceleration, hence, the input voltage of the laser vibrometer is not taken 
            into account. The default is 1.

        Returns
        -------
        None.

        """
        physics = laser_v_uff(file_path, rotation, self.frequency, useAverage, inputVoltage=inputVoltage)
        physics.ref2bem = ref2bem
        self.vibrometry[name] = physics
        self.radiator_id[name] = 'PLV'
        return None
    
    def vibrometry_data_user(self, name, Hv, X, ref2bem=None, inputVoltage=1):
        """
        Add acceleration data to the study. Meant to be used as a radiator - ref2bem strongly recommended.

        Parameters
        ----------
        name : str
            Reference.
        file_path : str
            Path to *.uff vibrometry data.
        rotation : list of float, optional
            Rotate measured data to match mesh direction. The default is [0, 0, 0].
        ref2bem : int, optional
            Reference to mesh physical group. The default is None.
        useAverage : Bool, optional
            Uses the average instead of individual element acceleration. 
            The default is False.
        inputVoltage : float, optional
            Optional scaling. By default, electroacPy uses the transfer function of the 
            acceleration, hence, the input voltage of the laser vibrometer is not taken 
            into account. The default is 1.

        Returns
        -------
        None.

        """
        physics = laser_v(Hv, X, inputVoltage)
        physics.ref2bem = ref2bem
        self.vibrometry[name] = physics
        self.radiator_id[name] = 'PLV'
        return None

    ## ========================
    # %% FILTERING / CROSSOVERS
    def filter_network(self, name, ref2bem=None, ref2study=None):
        """
        Add a filter network. If linked to a specific study, the correspondind surface should be linked in bem as well.
        :param name:
        :param ref2bem:
        :param ref2study:
        :return:
        """
        self.crossover[name] = xovers(self.frequency)
        self.crossover[name].ref2bem = ref2bem
        self.crossover[name].referenceStudy = ref2study
        return None

    def filter_addLowPass(self, networkRef, filterName, order, fc):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addLowPass(filterName, order, fc)
        else:
            self.crossover[networkRef].addLowPass(filterName, order, fc)
        return None

    def filter_addHighPass(self, networkRef, filterName, order, fc):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addHighPass(filterName, order, fc)
        else:
            self.crossover[networkRef].addHighPass(filterName, order, fc)
        return None

    def filter_addBandPass(self, networkRef, filterName, order, fc1, fc2):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addBandPass(filterName, order, fc1, fc2)
        else:
            self.crossover[networkRef].addBandPass(filterName, order, fc1, fc2)
        return None

    def filter_addLowPassBQ(self, networkRef, filterName, fc, Q, dBGain=0):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addLowPassBQ(filterName, fc, Q, dBGain)
        else:
            self.crossover[networkRef].addLowPassBQ(filterName, fc, Q, dBGain)
        return None

    def filter_addHighPassBQ(self, networkRef, filterName, fc, Q, dBGain=0):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addHighPassBQ(filterName, fc, Q, dBGain)
        else:
            self.crossover[networkRef].addHighPassBQ(filterName, fc, Q, dBGain)
        return None

    def filter_addBandPassBQ(self, networkRef, filterName, fc, Q, dBGain=0):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addBandPassBQ(filterName, fc, Q, dBGain)
        else:
            self.crossover[networkRef].addBandPassBQ(filterName, fc, Q, dBGain)
        return None

    def filter_addPeakEQ(self, networkRef, filterName, fc, Q, dBGain=0):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addPeakEQ(filterName, fc, Q, dBGain)
        else:
            self.crossover[networkRef].addPeakEQ(filterName, fc, Q, dBGain)
        return None

    def filter_addLowShelf(self, networkRef, filterName, fc, Q, dBGain=0):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addLowShelf(filterName, fc, Q, dBGain)
        else:
            self.crossover[networkRef].addLowShelf(filterName, fc, Q, dBGain)
        return None

    def filter_addHighShelf(self, networkRef, filterName, fc, Q, dBGain=0):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addHighShelf(filterName, fc, Q, dBGain)
        else:
            self.crossover[networkRef].addHighShelf(filterName, fc, Q, dBGain)
        return None

    def filter_addTransferFunction(self, networkRef, filterName, tf):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addTransferFunction(filterName, tf)
        else:
            self.crossover[networkRef].addTransferFunction(filterName, tf)
        return None

    def filter_addGain(self, networkRef, filterName, dBGain):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addGain(filterName, dBGain)
        else:
            self.crossover[networkRef].addGain(filterName, dBGain)
        return None

    def filter_addDelay(self, networkRef, filterName, dt, timeConvention='-jwt'):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addDelay(filterName, dt, timeConvention)
        else:
            self.crossover[networkRef].addDelay(filterName, dt, timeConvention)
        return None

    def filter_addPhaseFlip(self, networkRef, filterName):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].addPhaseFlip(filterName)
        else:
            self.crossover[networkRef].addPhaseFlip(filterName)
        return None

    def filter_delete(self, networkRef, filterName):
        if isinstance(networkRef, list):
            for net in networkRef:
                self.crossover[net].deleteFilter(filterName)
        else:
            self.crossover[networkRef].deleteFilter(filterName)
        return None

    ## ==================
    # %% ACOUSTIC STUDIES
    def study_acousticBEM(self, name, meshPath, acoustic_radiator, 
                          domain="exterior", **kwargs):
        """
        Add a boundary element analysis to the loudspeaker study.

        :param name:
        :param meshPath:
        :param acoustic_radiator: obj
            can either be an enclosure or driver object
        :param boundary:
        :param offset:
        :return:
        """
        rad_surf          = []
        surf_v            = []
        vibrometry_points = []
        Nfft              = len(self.frequency)
        if isinstance(acoustic_radiator, str):
            acoustic_radiator = [acoustic_radiator]   # set radiator in list is single str given

        # associate acoustic_radiator to velocity data // vibrometric data
        for i in range(len(acoustic_radiator)):
            cname = acoustic_radiator[i]

            # ========
            # enclosure
            if self.radiator_id[cname] == 'SPKBOX':
                try:
                    nRef = len(self.enclosure[cname].ref2bem)
                    for i in range(nRef):
                        rad_surf.append(self.enclosure[cname].ref2bem[i])
                        surf_v.append(np.ones(Nfft))
                        vibrometry_points.append(False)
                except:
                    rad_surf.append(self.enclosure[cname].ref2bem)
                    surf_v.append(np.ones(Nfft))
                    vibrometry_points.append(False)

            # ========
            # driver
            elif self.radiator_id[cname] == 'EAD':
                try:
                    nRef = len(self.driver[cname].ref2bem)
                    for i in range(nRef):
                        rad_surf.append(self.driver[cname].ref2bem[i])
                        surf_v.append(np.ones(Nfft))
                        vibrometry_points.append(False)
                except:
                    rad_surf.append(self.driver[cname].ref2bem)
                    surf_v.append(np.ones(Nfft))
                    vibrometry_points.append(False)

            # =========================
            # Polytech vibrometric data
            elif self.radiator_id[cname] == 'PLV':
                try:
                    nRef = len(self.vibrometry[cname].ref2bem)
                    for i in range(nRef):
                        rad_surf.append(self.vibrometry[cname].ref2bem[i])
                        surf_v.append(self.vibrometry[cname].v_point)
                        vibrometry_points.append(self.vibrometry[cname].point_cloud)
                except:
                    rad_surf.append(self.vibrometry[cname].ref2bem)
                    surf_v.append(self.vibrometry[cname].v_point)
                    vibrometry_points.append(self.vibrometry[cname].point_cloud)


        # surf_v = np.ones([len(rad_surf), len(self.frequency)])
        physics = bem(meshPath, rad_surf, surf_v, self.frequency, domain,
                      c_0=self.c, rho_0=self.rho, 
                      vibrometry_points=vibrometry_points, **kwargs)
        physics.radiator = acoustic_radiator
        self.acoustic_study[name] = physics
        self.evaluation[name] = evs_bem(physics)
        self.evaluation[name].referenceStudy = name
        return None

    
    def study_acousticPointSource(self, name, xSource, acoustic_radiator, 
                                  meshPath=None, domain="exterior", **kwargs):
        
        if isinstance(acoustic_radiator, str):
            acoustic_radiator = [acoustic_radiator]
        
        rad_point = []
        for cname in acoustic_radiator:
            if self.radiator_id[cname] == 'SPKBOX':
                try:
                    nRef = len(self.enclosure[cname].ref2bem)
                    for i in range(nRef):
                        rad_point.append(self.enclosure[cname].ref2bem[i])
                except:
                    rad_point.append(self.enclosure[cname].ref2bem)

            # ========
            # driver
            elif self.radiator_id[cname] == 'EAD':
                try:
                    nRef = len(self.driver[cname].ref2bem)
                    for i in range(nRef):
                        rad_point.append(self.driver[cname].ref2bem[i])
                except:
                    rad_point.append(self.driver[cname].ref2bem)
            
            
        if meshPath is not None:
            physics = pointSourceBEM(meshPath, xSource, 
                                     np.ones([len(xSource), 
                                              len(self.frequency)]), 
                                     self.frequency, domain=domain, 
                                     c_0=self.c, rho_0=self.rho, 
                                     radiatingElement=rad_point, 
                                     **kwargs)
            physics.radiator = acoustic_radiator
            self.acoustic_study[name] = physics
            self.evaluation[name] = evs_bem(physics)
            self.evaluation[name].referenceStudy = name
        else:
            physics = pointSource(xSource, np.ones([len(xSource), 
                                                    len(self.frequency)]), 
                                  self.frequency, radiatingElement=rad_point,
                                  c=self.c, rho=self.rho, **kwargs)
            physics.radiator = acoustic_radiator
            self.acoustic_study[name] = physics
            self.evaluation[name] = evs_bem(physics)
            self.evaluation[name].referenceStudy = name
        return None


    ## =======================
    # %% ACOUSTIC obs
    def evaluation_polarRadiation(self,
                           reference_study: str or list,
                           evaluation_name: str,
                           min_angle: float,
                           max_angle: float,
                           step: float,
                           on_axis: str,
                           direction: str,
                           radius: float = 5,
                           offset: list = [0, 0, 0]):
        """
        Add a circular microphone array to given study.

        :param reference_study:
        :param evaluation_name:
        :return:
        """
        if isinstance(reference_study, list):
            for i in range(len(reference_study)):
                self.evaluation[reference_study[i]].polarRadiation(evaluation_name,
                                                                    min_angle,
                                                                    max_angle,
                                                                    step,
                                                                    on_axis,
                                                                    direction,
                                                                    radius,
                                                                    offset)
        else:
            self.evaluation[reference_study].polarRadiation(evaluation_name,
                                                                min_angle,
                                                                max_angle,
                                                                step,
                                                                on_axis,
                                                                direction,
                                                                radius,
                                                                offset)
        return None

    def evaluation_pressureField(self,
                          reference_study: str,
                          evaluation_name: str,
                          L1: float,
                          L2: float,
                          step: float,
                          plane: str,
                          offset: list = [0, 0, 0],
                          addToPlotter: bool = True):
        """
        Add planar microphone array for "slice" pressure plot in given study.

        :param reference_study:
        :param evaluation_name:
        :param L1:
        :param L2:
        :param step:
        :param plane:
        :param offset:
        :param addToPlotter:
        :return:
        """

        if isinstance(reference_study, list):
            for i in range(len(reference_study)):
                self.evaluation[reference_study[i]].pressureField(evaluation_name,
                                                                   L1,
                                                                   L2,
                                                                   step,
                                                                   plane,
                                                                   offset)
        else:
            self.evaluation[reference_study].pressureField(evaluation_name,
                                                               L1,
                                                               L2,
                                                               step,
                                                               plane,
                                                               offset)
        return None

    def evaluation_fieldPoint(self,
                             reference_study: str,
                             evaluation_name: str,
                             microphonePositions: list,
                             **kwargs):
        """
        Add FRF evaluation points at defined position on given study.

        :param reference_study:
        :param evaluation_name:
        :param microphonePositions:
        :param labels:
        :return:
        """

        if isinstance(reference_study, list):
            for i in range(len(reference_study)):
                self.evaluation[reference_study[i]].fieldPoint(evaluation_name,
                                                                microphonePositions,
                                                                **kwargs)
        else:
            self.evaluation[reference_study].fieldPoint(evaluation_name,
                                                         microphonePositions,
                                                         **kwargs)
        return None

    def evaluation_boundingBox(self,
                                reference_study: str,
                                evaluation_name: str,
                                Lx: float,
                                Ly: float,
                                Lz: float,
                                step: float = 343/500/5,
                                offset: list = [0, 0, 0]):

        if isinstance(reference_study, list):
            for i in range(len(reference_study)):
                self.evaluation[reference_study[i]].boundingBox(evaluation_name,
                                                                 Lx, Ly, Lz, step, offset)
        else:
            self.evaluation[reference_study].boundingBox(evaluation_name,
                                                          Lx, Ly, Lz, step, offset)
        return None

    def evaluation_sphericalRadiation(self,
                                       reference_study: str,
                                       evaluation_name: str,
                                       nMic: float,
                                       radius: float = 1.8,
                                       offset: list = [0, 0, 0]):

        if isinstance(reference_study, list):
            for i in range(len(reference_study)):
                self.evaluation[reference_study[i]].sphericalRadiation(evaluation_name,
                                                                 nMic, radius, offset)
        else:
            self.evaluation[reference_study].sphericalRadiation(evaluation_name,
                                                          nMic, radius, offset)
        return None

    
    def evaluation_plottingGrid(self, 
                                reference_study: str, 
                                evaluation_name: str,
                                path_to_grid: str):
        
        if isinstance(reference_study, list):
            for i in range(len(reference_study)):
                self.evaluation[reference_study[i]].plottingGrid(evaluation_name,
                                                                 path_to_grid)
        else:
            self.evaluation[reference_study].plottingGrid(evaluation_name,
                                                          path_to_grid)

    ## ===================
    # %% run / plot / info
    def run(self):
        """
        Run all defined studies

        """

        for study in self.acoustic_study:
            if self.acoustic_study[study].isComputed is False:
                if hasattr(self.acoustic_study[study], "xSource") and self.acoustic_study[study].meshPath is None:
                    pass # check if it is a point source study without BEM boundaries. No "solve" if True
                else:
                    self.acoustic_study[study].solve()
            else:
                None
            
        for obs in self.evaluation:
            if bool(self.evaluation[obs].setup) is True:
                self.evaluation[obs].solve()
            else:
                print("No evaluation to compute for study {}.".format(obs))
        return None

    ## PLOT
    def plot_results(self, study=[], 
                     evaluation=[], radiatingElement=[], bypass_xover=False,
                     transformation="SPL", export_grid=False, pf2grid=False):
        
        # update solutions
        if isinstance(radiatingElement, int) is True: # avoid possible error if only one rad surf is selected
            radiatingElement = [radiatingElement]
                
        if bool(study) is False: # if not specific study given, plot all
            for s in self.acoustic_study:
                _ = updateResults(self, s, bypass_xover)
                _ = self.evaluation[s].plot(evaluation, radiatingElement,  
                                             processing=self.results[s],
                                             transformation=transformation,
                                             export_grid=export_grid, 
                                             pf2grid=pf2grid)
        else: # plot specific study
            _ = updateResults(self, study, bypass_xover)
            _ = self.evaluation[study].plot(evaluation, radiatingElement,
                                                 processing=self.results[study],
                                                 transformation=transformation,
                                                 export_grid=export_grid,
                                                 pf2grid=pf2grid)
        return None

    def plot_system(self, study, backend="pyvista"):
        """
        Plot study's mesh and related evaluations. By default it uses PyVista as the 
        plotting backend, but can use Gmsh if "gmsh" is passed.

        Parameters
        ----------
        study : str
            Study to display.
        backend : str, optional
            Which backend should be used. The default is "pyvista".

        Returns
        -------
        A plot with system mesh and evaluation points.

        """
        self.evaluation[study].plot_system(backend=backend)
        return None

    def plot_xovers(self, networks, split=True, amplitude=64):
        h_p = []
        names = []
        h = np.zeros(len(self.frequency), dtype=complex)
        if isinstance(networks, str) is True:
            self.crossover[networks].plot(divide_stages=split, amplitude=amplitude)
        else:
            for i in range(len(networks)):
                self.crossover[networks[i]].plot(divide_stages=split, amplitude=amplitude)
                h += self.crossover[networks[i]].h
                h_p.append(self.crossover[networks[i]].h)
                names.append(networks[i])
            plot_network(self, h, h_p, names, split, amplitude)
        return
    
    def plot_pressureMesh(self, study, radiatingElement=[],
                          bypass_xover=False, transformation="SPL", 
                          export_grid=False, backend="gmsh"):
        
        # update solutions
        if isinstance(radiatingElement, int) is True: # avoid possible error if only one rad surf is selected
            radiatingElement = [radiatingElement]
            
        if bool(radiatingElement) is False:
            # plot all elements
            element2plot = self.acoustic_study[study].radiatingElement
        else:
            element2plot = radiatingElement
                
        _ = updateResults(self, study, bypass_xover)        
        if backend == "gmsh":
            elementCoeff = np.ones([len(self.frequency), len(element2plot)], 
                                   dtype=complex)
            pp = self.results[study]
            for name in pp.TF:
                for idx, element in enumerate(element2plot):
                    if element in pp.TF[name]["radiatingElement"]:
                        elementCoeff[:, idx] *= pp.TF[name]["H"]
            
            bempp_grid_mesh(self.acoustic_study[study], elementCoeff, 
                            element2plot, transformation, export_grid)          
        return None


    ## Get Values
    # INFO
    def info(self):
        """
        Returns informations about the current simulation.

        Returns
        -------
        None
        """
        out_message = ()
        maxInfo = []
        if bool(self.driver) is True:
            ead_INFO = list(self.driver)
            out_message += (ead_INFO, )
            maxInfo.append(len(str(ead_INFO)))
        else:
            out_message += ('No driver defined', )
        if bool(self.enclosure) is True:
            box_INFO = list(self.enclosure)
            out_message += (box_INFO,)
            maxInfo.append(len(str(box_INFO)))
        else:
            out_message += ('No enclosure defined', )
        if bool(self.vibrometry) is True:
            laser_INFO = list(self.vibrometry)
            out_message += (laser_INFO, )
            maxInfo.append(len(str(laser_INFO)))
        else:
            out_message += ('No vibration data imported', )
        if bool(self.crossover) is True:
            xover_INFO = list(self.crossover)
            out_message += (xover_INFO,)
            maxInfo.append(len(str(xover_INFO)))
        else:
            out_message += ('No crossover defined',)
        if bool(self.acoustic_study) is True:
            study_INFO = list(self.acoustic_study)
            out_message += (study_INFO,)
            maxInfo.append(len(str(study_INFO)))
        else:
            out_message += ('No study defined',)
        
        maxInfo = np.max(maxInfo)
        maxDesc = len("# Electro-acoustic Drivers")
        maxTot = np.max([maxInfo,  maxDesc])

        print(maxTot * '#')
        print("# Electro-acoustic Drivers")
        print("# " + str(out_message[0]))

        print(maxTot * "-")
        print("# Enclosures")
        print("# " + str(out_message[1]))
        
        print(maxTot * '-')
        print("# Vibrometric data")
        print("# " + str(out_message[2]))

        print(maxTot * "-")
        print("# Crossovers")
        print("# " + str(out_message[3]))

        print(maxTot * "-")
        print("# Acoustical Studies")
        print("# " + str(out_message[4]))
        print(maxTot * "#")
        return None

    def get_pMic(self, studyName, evaluationName, 
                 radiatingElement=[], bypass_xover=False):
        """
        Return pressure at microphone points for defined studies, evaluation and radiating surfaces.

        :param studyName:
        :param evaluationName:
        :param radiatingSurface:
        :param get_freq_array:
        :return:
        """
        if bool(radiatingElement) is False:
            radiatingElement = self.acoustic_study[studyName].radiatingElement 
        
        if isinstance(radiatingElement, int) is True: # avoid possible error if only one rad surf is selected
            radiatingElement = [radiatingElement]

        _  = updateResults(self, studyName, bypass_xover)
        pmic = self.evaluation[studyName].setup[evaluationName].pMic
        
        if self.results[studyName] is not None:
            elementCoeff = np.ones([len(self.frequency), len(radiatingElement)], 
                                   dtype=complex)
            pp = self.results[studyName]
            for name in pp.TF:
                if bypass_xover is False:
                    for idx, element in enumerate(radiatingElement):
                        if element in pp.TF[name]["radiatingElement"]:
                            elementCoeff[:, idx] *= pp.TF[name]["H"]
                elif bypass_xover is True and name[:12] != "filter_stage":
                    for idx, element in enumerate(radiatingElement):
                        if element in pp.TF[name]["radiatingElement"]:
                            elementCoeff[:, idx] *= pp.TF[name]["H"]
        else:
            elementCoeff = np.ones([len(self.frequency), len(radiatingElement)], 
                                   dtype=complex)
        
        out = getPressure(pmic, self.acoustic_study[studyName].radiatingElement,
                                radiatingElement, elementCoeff)
        return out
    
    
    def export_directivity(self, folder_path, file_name,
                           study, evaluation, radiatingElement=[], 
                           bypass_xover=False, frd=False):
        """
        Export directivity results.

        Parameters
        ----------
        folder_name : str
            where data are saved.
        file_name: str
            file prefix.
        study : str
            name of study to export.
        evaluation : str
            evaluation to export.
        radiatingElement : int or list of int, optional
            extract specific element.
        bypass_xover : bool, optional 
            if True, will export filtered response.
        frd : bool, optional
            if True, uses *.frd file extension instead of *.txt

        Returns
        -------
        None.

        """
        from electroacPy.general.acoustics import export_directivity
        
        pmic = self.get_pMic(study, evaluation, radiatingElement, bypass_xover)
        theta = self.evaluation[study].setup[evaluation].theta
        frequency = self.frequency
        export_directivity(folder_path, file_name, frequency, theta, pmic, frd)
        
        
    def export_impedance(self, folder_path, file_name, objName, zma=False):
        """
        Export impedance into .txt file

        Parameters
        ----------
        folder_name : str
            where data are saved.
        file_name: str
            file prefix.
        objName : str
            enclosure or driver object to export.
        zma : bool, optional
            if True, uses the *.zma file extension.
            
        Returns
        -------
        None.

        """
        
        if zma is False:
            if objName in self.enclosure:
                self.enclosure[objName].exportZe(folder_path, file_name + ".txt")
            elif objName in self.driver:
                self.driver[objName].exportZe(folder_path, file_name + ".txt")
        elif zma is True:
            if objName in self.enclosure:
                self.enclosure[objName].exportZe(folder_path, file_name + ".zma")
            elif objName in self.driver:
                self.driver[objName].exportZe(folder_path, file_name + ".zma")
        
    
## =======================
# %% post-processing tools
def updateResults(loudspeakerSystem, study_to_plot, bypass_xover):
    ls = loudspeakerSystem
    study = study_to_plot
    ls.results[study] = pp()
    radiatorName      = ls.acoustic_study[study].radiator
    for name in radiatorName:
        if ls.radiator_id[name] == 'SPKBOX':
            _ = apply_Velocity_From_SPKBOX(ls, study, name)
        elif ls.radiator_id[name] == 'EAD':
            _ = apply_Velocity_From_EAD(ls, study, name)
        elif ls.radiator_id[name] == 'PLV':
            _ = apply_Velocity_From_PLV(ls, study, name)

    if bypass_xover is False:
        for xover in ls.crossover:
            if study in ls.crossover[xover].referenceStudy:  # check if some of the crossovers are to be applied in the study
                h = ls.crossover[xover].h
                ref2bem = ls.crossover[xover].ref2bem
                if isinstance(ref2bem, int):
                    ls.results[study].addTransferFunction("filter_stage_"+xover, 
                                                          h, [ref2bem])
                elif isinstance(ref2bem, list):
                    ls.results[study].addTransferFunction("filter_stage_"+xover, 
                                                          h, ref2bem)
                else:
                    print("you shouldn't use a numpy array, at least for now :)")
    return None

def apply_Velocity_From_SPKBOX(loudspeakerSystem, study, radiatorName):
    """
    Apply velocity from speakerBox object on BEM simulation results

    :param loudspeakerSystem:
    :param radiatorName:
    :param portOnly
    :return:
    """
    ls = loudspeakerSystem
    if "xSource" not in ls.acoustic_study[study].__dict__: # check if point_source
        v = ls.enclosure[radiatorName].v      # driver velocity
        vp = ls.enclosure[radiatorName].vp    # port's velocity
        vp2 = ls.enclosure[radiatorName].vp2  # port 2 velocity (in case of bp2 enclosure config)
        vpr = ls.enclosure[radiatorName].vpr  # passive radiator velocity (in case of pr enclosure config)
        vpr2 = ls.enclosure[radiatorName].vpr2
    else:
        v = ls.enclosure[radiatorName].Q     # driver velocity
        vp = ls.enclosure[radiatorName].Qp    # port's velocity
        vp2 = ls.enclosure[radiatorName].Qp2  # port 2 velocity (in case of bp2 enclosure config)
        vpr = ls.enclosure[radiatorName].Qpr  # passive radiator velocity (in case of pr enclosure config)
        vpr2 = ls.enclosure[radiatorName].Qpr2
    ref2bem = ls.enclosure[radiatorName].ref2bem
        
    if isinstance(ref2bem, int) is True:
        if ls.enclosure[radiatorName].config == "sealed":  # check if 1 radiator because sealed or bp config
            ls.results[study].addTransferFunction("driver_"+radiatorName,
                                                  v, [ref2bem])
        elif ls.enclosure[radiatorName].config == "bandpass":
            ls.results[study].addTransferFunction("port_"+radiatorName, 
                                                  vp, [ref2bem])
        elif ls.enclosure[radiatorName].config == "bandpass_pr":
            ls.results[study].addTransferFunction("passive-radiator_"+"radiatorName",
                                                  vpr, [ref2bem])
    elif isinstance(ref2bem, list) is True:
        if ls.enclosure[radiatorName].config == "sealed":
            ls.results[study].addTransferFunction("driver_"+radiatorName, v, ref2bem)
        elif ls.enclosure[radiatorName].config == "vented":
            ls.results[study].addTransferFunction("driver_"+radiatorName, v, ref2bem[:-1])
            ls.results[study].addTransferFunction("port_"+radiatorName, vp, [ref2bem[-1]])
        
        elif ls.enclosure[radiatorName].config == "passiveRadiator":
            ls.results[study].addTransferFunction("driver_"+radiatorName, v, ref2bem[:-1])
            ls.results[study].addTransferFunction("passive-radiator_"+radiatorName, vpr, [ref2bem[-1]])
        
        elif ls.enclosure[radiatorName].config == "bandpass":
            ls.results[study].addTransferFunction("port_"+radiatorName, vp, ref2bem)
        
        elif ls.enclosure[radiatorName].config == "bandpass_2":
            ls.results[study].addTransferFunction("portf_"+radiatorName, vp, [ref2bem[0]])
            ls.results[study].addTransferFunction("portb_"+radiatorName, vp2, [ref2bem[1]])

        elif ls.enclosure[radiatorName].config == "bandpass_pr":
            ls.results[study].addTransferFunction("passive-radiator_"+radiatorName, vpr, ref2bem)
        
        elif ls.enclosure[radiatorName].config == "bandpass_pr_2":
            ls.results[study].addTransferFunction("prf_"+radiatorName, vpr, [ref2bem[0]])
            ls.results[study].addTransferFunction("prb_"+radiatorName, vpr2, [ref2bem[1]])        
    return None

def apply_Velocity_From_EAD(loudspeakerSystem, study, radiatorName):
    """
    Apply velocity from EAD object on BEM simulation results

    :param loudspeakerSystem:
    :param radiatorName:
    :return:
    """
    ls = loudspeakerSystem
    
    if "xSource" not in ls.acoustic_study[study].__dict__:
        v = ls.driver[radiatorName].v
    else:
        v = ls.driver[radiatorName].Q
    ref2bem = ls.driver[radiatorName].ref2bem
    if isinstance(ref2bem, int) is True:
        ls.results[study].addTransferFunction("driver_"+radiatorName, v, [ref2bem])
    elif isinstance(ref2bem, list) is True:
        ls.results[study].addTransferFunction("driver_"+radiatorName, v, ref2bem)
    return None

def apply_Velocity_From_PLV(loudspeakerSystem, study, radiatorName):
    """
    Apply decoy velocity from PLV object on BEM simulation results

    :param loudspeakerSystem:
    :param radiatorName:
    :return:
    """
    ls = loudspeakerSystem
    v = ls.vibrometry[radiatorName].v
    ref2bem = ls.vibrometry[radiatorName].ref2bem
    if isinstance(ref2bem, int) is True:
        ls.results[study].addTransferFunction("unit_v_"+radiatorName, v, [ref2bem])
    elif isinstance(ref2bem, list) is True:
        ls.results[study].addTransferFunction("unit_v_"+radiatorName, v, ref2bem)
    return None

## ===============
# %% plot function
def plot_network(system, h, h_p, names, split, amplitude):
    if split is False:
        fig, ax = plt.subplots(2, 1)
        ax[0].semilogx(system.frequency, dB(h), 'k')
        ax[1].semilogx(system.frequency, np.angle(h), 'k')
        ax[0].set(xlabel='Frequency [Hz]', ylabel='Gain [dB]',
               title='Total network response',
               ylim=[-amplitude/2, +amplitude/2])
        ax[1].set(xlabel='Frequency [Hz]', ylabel='Phase [rad]')
        for i in range(2):
            ax[i].grid(linestyle='dotted', which='both')
            ax[i].legend(loc='best')
        plt.tight_layout()
    else:
        fig, ax = plt.subplots(2, 1)
        for i in range(len(h_p)):
            ax[0].semilogx(system.frequency, dB(h_p[i]), label=names[i])
            ax[1].semilogx(system.frequency, np.angle(h_p[i]), label=names[i])
        ax[0].semilogx(system.frequency, dB(h), 'k', label='summed response')
        ax[1].semilogx(system.frequency, np.angle(h), 'k', label='summed response')
        ax[0].set(xlabel='Frequency [Hz]', ylabel='Gain [dB]',
               title='Total network response',
               ylim=[-amplitude/2, +amplitude/2])
        ax[1].set(xlabel='Frequency [Hz]', ylabel='Phase [rad]')
        for j in range(2):
            ax[j].grid(linestyle='dotted', which='both')
            ax[j].legend(loc='best')
        plt.tight_layout()
    return plt.show()
