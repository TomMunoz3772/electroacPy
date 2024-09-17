#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:29:39 2024

@author: tom
"""

import numpy as np
import generalToolbox as gtb
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Slider, TextBox, Button
from matplotlib.backend_bases import MouseEvent, MouseButton
from matplotlib.contour import QuadContourSet
from matplotlib import cm, colors
import vtk
import pyvista
import os
import electroacPy
directory_path = os.path.abspath(electroacPy.__file__)

pi = np.pi

class observations:
    """
    Tools for creating and managing various types of acoustic observations using a Boundary Element Method (BEM) object.

    Available observation types:
    - Planar observations: Defined on a planar surface with specified dimensions.
    - Polar observations: Defined along a circular path in a specified plane.
    - Spherical observations: Defined on a spherical surface.
    - Pressure response observations: Recorded pressure responses at specified microphone positions.

    Parameters:
        bemObject (BEM object): A Boundary Element Method (BEM) object after computing pressure on the mesh.

    Attributes:
        identifier (str): Class identifier.
        bemObject: The BEM object used for calculations.
        observationName (list): Names of added observations.
        observationType (list): Types of added observations.
        computedObservations (list): Names of computed observations.
        xMic (list): Microphone positions for each observation.
        labels (list): Labels associated with observations.
        L (list): Planar observation lengths.
        W (list): Planar observation widths.
        planarPlane (list): Planes on which planar observations are defined.
        theta (list): Angles for polar observations.
        polarPlane (list): Planes on which polar observations are defined.
        DI (list): Directivity indices of polar observations.
        pMic (list): Pressure data for each computed observation.
        pMicArray (list): Pressure data for each computed observation with separate speakers.

    Methods:
        addPlanarObservation: Add a planar observation setup.
        addPolarObservation: Add a polar observation setup.
        addPressureResponse: Add a pressure response observation setup.
        addSphericalObservation: Add a spherical observation setup.
        computeObservations: Compute the added observations.
        plot: Plot the computed observations.
        plot_system: Plot the BEM mesh and observation points in 3D space.
        get_pMic: Retrieve the pressure data of a specific observation.
        get_directivityIndex: Compute and plot the directivity index of a polar observation.
        deleteObservation: Delete a specified observation.
        save: Save computed observations and related data to an .npz file.
    """
    def __init__(self, bemObject):
        self.bemObject = bemObject
        self.setup = {}
        
        # ref to system
        self.referenceStudy = None
    
    
    def polarRadiation(self, observationName, minAngle: float, maxAngle: float,
                       step: float, on_axis: str, direction: str, 
                       radius: float = 5, offset: list = [0, 0, 0], **kwargs):
        
        if observationName not in self.setup:
            self.setup[observationName] = PolarRadiation(minAngle, maxAngle, 
                                                         step, on_axis, direction,
                                                         offset)
        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = PolarRadiation(minAngle, maxAngle, 
                                                         step, on_axis, direction,
                                                         offset)
        else:
             print("observation {} already exists. You can overwrite it using  \
                   the 'overwrite' flag".format(observationName))   
            
    
    def pressureField(self, observationName, Length: float, Width: float,
                      step: float, plane: str, offset: list = [0, 0, 0],
                      **kwargs):
        
        if observationName not in self.setup:
            self.setup[observationName] = PressureField(Length, Width, step, 
                                                        plane, offset)
        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = PressureField(Length, Width, step, 
                                                        plane, offset)
        else:
             print("observation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(observationName))

        
    def fieldPoint(self, observationName, microphonePosition, **kwargs):        
        if observationName not in self.setup:
            self.setup[observationName] = FieldPoint(microphonePosition, kwargs)

        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = FieldPoint(microphonePosition, kwargs)

        else:
             print("observation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(observationName))   


    def boundingBox(self, observationName, Lx, Ly, Lz,
                    step=1, offset=[0, 0, 0], **kwargs):
        if observationName not in self.setup:
            self.setup[observationName] = BoundingBox(Lx, Ly, Lz, step, offset)

        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = BoundingBox(Lx, Ly, Lz, step, offset)

        else:
             print("observation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(observationName))   
        
        
    def sphericalRadiation(self, observationName, nMic, 
                           radius, offset=[0, 0, 0], **kwargs):
        if observationName not in self.setup:
            self.setup[observationName] = SphericalRadiation(nMic, radius, offset)

        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = SphericalRadiation(nMic, radius, offset)

        else:
             print("observation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(observationName))  
    
    def computeObservations(self, observation_name="all"):
        if observation_name == "all":
            obs_to_compute = list(self.setup.keys())
        elif isinstance(observation_name, str):
            obs_to_compute = [observation_name]
        elif isinstance(observation_name, list):
            obs_to_compute = observation_name
        else:
            raise ValueError("observation_name should be a key in setup.")
        
        # group microphones to compute in one step (faster for large numbers
        # of microphones)
        mic2compute = np.empty([0, 3])
        for obs in obs_to_compute:
            setup = self.setup[obs]
            xMic = setup.xMic
            nMic = setup.nMic
            if setup.isComputed is False:
                mic2compute = np.concatenate((mic2compute, xMic))
        
        _, pMic = self.bemObject.getMicPressure(mic2compute, 
                                                individualSpeakers=True)
        
        # ungroup microphones and store within each setup
        current_index = 0
        for obs in obs_to_compute:
            setup = self.setup[obs]
            nMic = self.nMic
            setup.pMic = pMic[:, current_index:current_index+nMic, :]
            setup.isComputed = True

            
            
            
            
            
            
    
#%% Evaluation classes
class EvaluationParameters:
    def __init__(self):
        self.type = None
        self.xMic = None
        self.isComputed = False
        self.pMic = None
        
        
class PolarRadiation(EvaluationParameters):
    def __init__(self, minAngle: float,
                 maxAngle: float,
                 step: float,
                 on_axis: str,
                 direction: str,
                 radius: float = 5,
                 offset: list = [0, 0, 0]):
        self.minAngle = minAngle
        self.maxAngle = maxAngle
        self.step = step
        self.on_axis = on_axis
        self.direction = direction
        self.radius = radius
        self.offset = offset
        
        from generalToolbox.geometry import create_circular_array
        self.theta = np.arange(minAngle, maxAngle+step, step)
        self.xMic = create_circular_array(self.theta, 
                                          on_axis, direction, radius, offset)
        self.nMic = len(self.xMic)
        # edit evaluationParameters
        self.type = "polar"
        

class FieldPoint(EvaluationParameters):
    def __init__(self, micPositions, **kwargs): 
        self.xMic = np.array(micPositions)
        self.nMic = len(self.xMic)

        if "legend" in kwargs:
            self.legend = kwargs["legend"]
        
        # edit evaluationParameters
        self.type = "field_point"
        
        
class PressureField(EvaluationParameters):
    def __init__(self, Length: float, 
                 Width: float,
                 step: float,
                 plane: str,
                 offset: list = [0, 0, 0]):
        self.Width = Width
        self.Length = Length
        self.step = step
        self.plane = plane
        self.offset = offset
        
        from generalToolbox.geometry import create_planar_array
        self.xMic, self.L, self.W = create_planar_array(Length, Width, step, plane)        
        self.nMic = len(self.xMic)

        # edit evaluationParameters
        self.type = "pressure_field"

class BoundingBox(EvaluationParameters):
    def __init__(self, Lx, Ly, Lz, step=1, offset=[0, 0, 0]):
        self.Lx = Lx
        self.Ly = Ly 
        self.Lz = Lz
        self.step = step
        self.offset = offset
        
        from generalToolbox.geometry import create_bounding_box
        self.xMic = create_bounding_box(Lx, Ly, Lz, step, offset)
        self.nMic = len(self.xMic)

        # edit evaluationParameters
        self.type = "bounding_box"
        
        
class SphericalRadiation(EvaluationParameters):
    def __init__(self, nMic, radius, offset):
        self.radius = radius
        self.offset = offset
        
        from generalToolbox.geometry import create_spherical_array
        self.xMic = create_spherical_array(nMic, radius, offset)
        self.nMic = len(self.xMic)

        # edit evaluationParameters
        self.type = "spherical"
        
        
        
        
        
        