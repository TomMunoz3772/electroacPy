#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:43:24 2023

@author: tom.munoz
"""

import numpy as np
from electroacPy.acousticSim.interiorAcoustics.observations import observations
from electroacPy.acousticSim.interiorAcoustics.bem import bem as bem_int
from copy import copy, deepcopy
import bempp.api

class postProcess:
    """
    Tool for post-processing acoustic observations by applying filters to radiating surfaces.

    This class is designed to work with an existing `observations` object, adding filtering capabilities
    to modify the pressure data associated with different radiating surfaces.

    Parameters:
        observationObject (observations object): The observation object containing computed pressure data.

    Attributes:
        identifier (str): Class identifier.
        observationObject: The `observations` object to which post-processing filters will be applied.
        freq_array: Frequency array from the associated `observations` object.
        radiatingSurface: List of radiating surfaces from the associated `observations` object.
        observationName: List of computed observation names.
        add2Plotter: List of flags indicating whether observations should be added to the plotter.
        pMicArray: Deep copy of the pressure data arrays for each observation and radiating surface.

    Methods:
        addFilter: Apply a transfer function filter to specific radiating surfaces.
        plot: Plot the post-processed observations.
        get_pMic: Retrieve the pressure data of a specific post-processed observation.
        save: Save the post-processed observations and related data to a file.
    """
    def __init__(self, observationObject):
        self.identifier = 'PostP'
        self.observationObject = observationObject

        # parameters
        self.freq_array = observationObject.bemObject.freq_array
        self.radiatingSurface = observationObject.bemObject.radiatingSurface
        self.observationName = observationObject.computedObservations
        self.add2Plotter = observationObject.add2Plotter

        # mic data to update (initialisation)
        self.pMicArray = deepcopy(observationObject.pMicArray)
        # self.pMeshArray = deepcopy(observationObject.bemObject.p_total_list_array[:][0].coefficients)

    def addFilter(self, transferFunction, radiatingSurface):
        """
        Apply the given transfer function to the specified radiating surface.

        Parameters:
        ----------
        transferFunction : ndarray
            Transfer function of complex filter. The frequency array on which it is defined must correspond to the observation's frequency array.
        radiatingSurface : list or ndarray
            Indices corresponding to the surface indices given in the mesh used for Boundary Element Method (BEM).

        Returns:
        -------
        None
            This method updates the internal state of the observation object to include the effect of the applied filter on the specified radiating surfaces.

        Notes:
        -----
        - The method applies the given transfer function to each radiating surface specified by their indices in the provided list or ndarray.
        - The observation's radiating surfaces are assumed to be pre-defined in the observation object.
        - The updated observation object can be accessed using the observationUpdate attribute after calling this method.
        """
        radSurf_bem = self.radiatingSurface # object radsurface
        pMesh       = np.empty([len(self.freq_array), len(self.radiatingSurface)], dtype=object)
        nRadSurf    = len(radSurf_bem)
        pMicArray_update =  self.pMicArray # copy of pMicArray
        pMic_update = []

        # microphone pressure
        indexSurface = []
        for pMicIndex in range(len(pMicArray_update)):
            for i, val in enumerate(radSurf_bem):
                if val in radiatingSurface:
                    dimA = np.shape(pMicArray_update[pMicIndex][:, :, i])[0]
                    dimB = np.shape(pMicArray_update[pMicIndex][:, :, i])[1]
                    TF = np.zeros([dimA, dimB], dtype=complex)
                    for k in range(dimB):
                        TF[:, k] = transferFunction
                    pMicArray_update[pMicIndex][:, :, i] *= TF
            pMic_update.append(np.sum(pMicArray_update[pMicIndex], 2))

        # mesh pressure
        for f in range(len(self.freq_array)):
            for i, val in enumerate(radSurf_bem):
                coeff = deepcopy(self.observationObject.bemObject.p_mesh[f, i].coefficients)
                if val in radiatingSurface:
                    coeff *= transferFunction[f]
                p_total_mesh = bempp.api.GridFunction(self.observationObject.bemObject.spaceP, coefficients=coeff)
                pMesh[f, i] = p_total_mesh

        # store updated microphone data (should be updated as more filters are added)
        self.pMicArray = pMicArray_update

        # set updated dataset in new observation object
        self.observationUpdate = buildPostProcessedObservation(self.observationObject,
                                                               pMicArray_update,
                                                               pMic_update,
                                                               pMesh)
        return None

    def plot(self, obs_name='all', radiatingSurface='all'):
        return self.observationUpdate.plot(obs_name, radiatingSurface)

    def get_pMic(self, obsName, radiatingSurface='all', get_freq_array=False):
        return self.observationUpdate.get_pMic(obsName, radiatingSurface, get_freq_array)

    def save(self, filename):
        return self.observationUpdate.save(filename)


def buildPostProcessedObservation(observationObject, # observation data (not processed)
                                  pMicArray, pMic, pMesh):  # filtered data
    """
    Create a new observation object from the updated pressure responses.

    Parameters:
    ----------
    observationObject : observations
        Original observation data (not processed).
    pMicArray : list of ndarray
        Updated pressure microphone array data after applying filters.
    pMic : list of ndarray
        Updated processed pressure microphone data after applying filters.

    Returns:
    -------
    obsUpdate : observations
        New observation object containing updated pressure responses.

    Notes:
    -----
    - This function creates a new observation object based on the provided original observation data and the updated pressure microphone array and processed data.
    - The new observation object inherits various properties from the original observationObject but with updated pressure data.
    - The other attributes like observationName, observationType, computedObservations, xMic, labels, planar properties, polar properties, etc., are copied from the original observationObject.
    """

    meshPath            = observationObject.bemObject.meshPath
    radSurf             = observationObject.bemObject.radiatingSurface
    freq_array          = observationObject.bemObject.freq_array
    surfVelocity        = observationObject.bemObject.surfaceVelocity
    vibrometry_points   = observationObject.bemObject.vibrometry_points
    radiation_direction = observationObject.bemObject.radiation_direction
    absSurface          = observationObject.bemObject.absorbingSurface
    surfImp             = observationObject.bemObject.surfaceImpedance
    c_0                 = observationObject.bemObject.c_0
    rho_0               = observationObject.bemObject.rho_0

    bemObject = bem_int(meshPath, radSurf, surfVelocity, freq_array, absSurface, surfImp,
                        vibrometry_points, radiation_direction, c_0=c_0, rho_0=rho_0)

    bemObject.p_mesh = pMesh
    obsUpdate = observations(bemObject)  # new object

    obsUpdate.observationName = copy(observationObject.observationName)
    obsUpdate.observationType = copy(observationObject.observationType)
    obsUpdate.computedObservations = copy(observationObject.computedObservations)

    # xmic
    obsUpdate.xMic = copy(observationObject.xMic)
    obsUpdate.labels = copy(observationObject.labels)

    # planar
    obsUpdate.L = copy(observationObject.L)
    obsUpdate.W = copy(observationObject.W)
    obsUpdate.planarPlane = copy(observationObject.planarPlane)
    obsUpdate.add2Plotter = copy(observationObject.add2Plotter)

    # polar
    obsUpdate.theta = copy(observationObject.theta)
    obsUpdate.polarPlane = copy(observationObject.polarPlane)
    obsUpdate.DI = []

    # bounding box
    # print(observationObject.Lx)
    # print(observationObject.Ly)
    # print(observationObject.Lz)
    # obsUpdate.Lx = copy(observationObject.Lx)
    # obsUpdate.Ly = copy(observationObject.Ly)
    # obsUpdate.Lz = copy(observationObject.Lz)

    # simulated pressure
    obsUpdate.pMicArray = pMicArray
    obsUpdate.pMic = pMic
    return obsUpdate
