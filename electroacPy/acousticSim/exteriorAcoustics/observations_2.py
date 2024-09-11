#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:43:16 2023

@author: tom.munoz
"""

# import numpy as np
# import generalToolbox as gtb
# import matplotlib.pyplot as plt
# from matplotlib.widgets import Cursor, Slider, TextBox, Button
# from matplotlib.backend_bases import MouseEvent, MouseButton
# from matplotlib.contour import QuadContourSet
# from matplotlib import cm, colors
# import vtk
# import pyvista
# # from pyvistaqt import BackgroundPlotter
# import os
# import electroacPy

import numpy as np
import matplotlib.pyplot as plt
from generalToolbox.geometry import (createCircArray, createPlaneArray,
                                     createSphereArray, create_boundingBox)
from generalToolbox.plot import plotTotalDirectivity

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

    def __init__(self):
        self.identifier = "OBS_BEM_EXT"

        # After computation
        self.isComputed = False
        self.pressure = None

    def compute(self):
        """
        Evaluate the pressure at specified observations. Will group all observation define as input and do a single run.
        Don't hesitate to set a lot of observations, the computation time does not increase linearly. Three cases can
        be run:
        - all observation (default)
        - list of observations
        -
        :param obs_name: list, str,
            name of observations to evaluate
        :return:
        """

        if self.isComputed is False:
            _, self.pressure = self.bemObject.getMicPressure(self.xMic, individualSpeakers=True)
            self.isComputed = True

    def plotSystem(self):
        return None

class polarRadiation(observations):
    def __init__(self, bemObject,
                       obsName: str,
                       minAngle: float,
                       maxAngle: float,
                       step: float,
                       plane: str,
                       radius: float = 5,
                       offset: list = [0, 0, 0]):
        super().__init__()

        # init object
        self.type = "polar"
        self.name = obsName
        self.bemObject = bemObject
        self.freq_array = bemObject.freq_array

        self.theta = np.deg2rad(np.arange(minAngle, maxAngle + step, step))
        self.xmic = createCircArray(theta, plane, radius=radius, offset=offset)

    def plotResult(self):
        if self.isComputed is True:
            plotWindow = plotTotalDirectivity(self.theta, self.freq_array, )

