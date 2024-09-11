#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:43:16 2023

@author: tom.munoz
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
# from pyvistaqt import BackgroundPlotter
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
        # identifier
        self.identifier = "OBS_BEM_EXT"

        self.bemObject = bemObject
        self.observationName = []
        self.observationType = []
        self.computedObservations = []

        # xmic
        self.xMic = []
        self.labels = []

        # planar
        self.L = []
        self.W = []
        self.planarPlane = []
        self.add2Plotter = []

        # polar
        self.theta = []
        self.polarPlane = []
        self.DI = []

        # bounding box -> not sure that it will be used
        self.Lx = []
        self.Ly = []
        self.Lz = []

        # simulated pressure
        self.pMic = []        # all speakers playing
        self.pMicArray = []   # pMic with separate speakers
        
        # ref to system
        self.referenceStudy = None
        
        
    def pressureField(self, obsName: str,
                             Length: float, 
                             Width: float,
                             step: float,
                             plane: str,
                             offset: list = [0, 0, 0],
                             addToPlotter=False):

        # check if observation exists
        isComputed = self.checkObservationExists(obsName, deleteObservation=True)

        if isComputed is False:
            # create mic array
            vert = self.bemObject.grid_show.vertices.T
            xmic, L, W = createPlaneArray(Length, Width, step, plane, offset)

            # update observations
            self.observationName.append(obsName)
            self.observationType.append('planar')

            # update observation arrays
            self.xMic.append(xmic)

            # planar
            self.L.append(L)
            self.W.append(W)
            self.planarPlane.append(plane)
            self.add2Plotter.append(addToPlotter)

            # polar -- 'trash' values
            self.theta.append('N/A')
            self.polarPlane.append('N/A')

            # frf
            self.labels.append('N/A')

            # bounding box
            self.Lx.append('N/A')
            self.Ly.append('N/A')
            self.Lz.append('N/A')
        return None
    
    
    def polarRadiation(self, obsName: str,
                             minAngle: float,
                             maxAngle: float,
                             step: float,
                             # plane: str,
                             on_axis: str,
                             direction: str,
                             radius: float = 5,
                             offset: list = [0, 0, 0]):

        # check if observation exists
        isComputed = self.checkObservationExists(obsName, deleteObservation=True)

        if isComputed is False:
            # create mic array
            from generalToolbox.geometry import create_circular_array_2

            theta = np.arange(minAngle, maxAngle+step, step)
            # xmic = createCircArray(theta, plane, radius=radius, offset=offset)
            xmic = create_circular_array_2(theta, on_axis, direction, radius, offset)

            # update observation names and types
            self.observationName.append(obsName)
            self.observationType.append('polar')

            # update observation arrays
            self.xMic.append(xmic)

            # polar
            self.theta.append(np.deg2rad(theta))
            # self.polarPlane.append(plane)
            self.polarPlane.append(on_axis+direction)

            # planar -- "trash" values
            self.L.append('N/A')
            self.W.append('N/A')
            self.planarPlane.append('N/A')
            self.add2Plotter.append(False)

            # frf
            self.labels.append('N/A')

            # bounding box
            self.Lx.append('N/A')
            self.Ly.append('N/A')
            self.Lz.append('N/A')
        return None
    
    
    def pressureResponse(self, obsName: str,
                            micPositions: list,
                            labels: list):

        # check if observation exists
        isComputed = self.checkObservationExists(obsName, deleteObservation=True)

        if isComputed is False:
            xmic = np.array(micPositions)

            self.xMic.append(xmic)
            self.observationName.append(obsName)
            self.observationType.append('FRF')

            # planar -- "trash" values
            self.L.append('N/A')
            self.W.append('N/A')
            self.planarPlane.append('N/A')
            self.add2Plotter.append(False)

            # polar -- 'trash' values
            self.theta.append('N/A')
            self.polarPlane.append('N/A')

            # frf
            self.labels.append(labels)

            # bounding box
            self.Lx.append('N/A')
            self.Ly.append('N/A')
            self.Lz.append('N/A')
        return None


    def boundingBox(self,  obsName: str,
                        Lx: float,
                        Ly: float,
                        Lz: float,
                        step: float = 343/500/5,
                        offset: list = [0, 0, 0]):

        isComputed = self.checkObservationExists(obsName, deleteObservation=True)

        if isComputed is False:
            xmic = create_boundingBox(Lx, Ly, Lz, step, offset)
            self.xMic.append(xmic)
            self.observationName.append(obsName)
            self.observationType.append('BOX')

            # planar -- "trash" values
            self.L.append('N/A')
            self.W.append('N/A')
            self.planarPlane.append('N/A')
            self.add2Plotter.append(False)

            # polar -- 'trash' values
            self.theta.append('N/A')
            self.polarPlane.append('N/A')

            # frf
            self.labels.append('N/A')

            # bounding box
            self.Lx.append(int((Lx+step) / step) + 1)
            self.Ly.append(int((Ly+step) / step) + 1)
            self.Lz.append(int((Lz+step) / step) + 1)
        return None


    def sphericalRadiation(self, obsName: str,
                           nMic: int = 256,
                           radius: float = 1.8,
                           offset: list = [0, 0, 0]):
        # check if observation exists
        isComputed = self.checkObservationExists(obsName, deleteObservation=True)

        if isComputed is False:
            xmic = createSphereArray(nMic, radius)

            # update observation names and types
            self.observationName.append(obsName)
            self.observationType.append('spherical')

            # update observation arrays
            self.xMic.append(xmic)

            # polar -- "trash" values
            self.theta.append('N/A')
            self.polarPlane.append('N/A')

            # planar -- "trash" values
            self.L.append('N/A')
            self.W.append('N/A')
            self.planarPlane.append('N/A')
            self.add2Plotter.append(False)

            # frf
            self.labels.append('N/A')

            # bounding box
            self.Lx.append('N/A')
            self.Ly.append('N/A')
            self.Lz.append('N/A')
        return None

    def computeObservations(self, obs_name='all'):
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
        micToCompute = np.empty([0, 3])
        micNumber = []
        if obs_name == 'all':
            nObs = len(self.observationName)
            for i in range(nObs):
                isComputed = self.checkObservationExists(self.observationName[i])
                if isComputed is False:
                    micToCompute = np.concatenate((micToCompute, self.xMic[i]), axis=0)
                    self.computedObservations.append(self.observationName[i])
                    micNumber.append(len(self.xMic[i]))
            print('\n' + '=' * 40)
            # print('Solving pressure for all observations')
            pMic, pMicArray = self.bemObject.getMicPressure(micToCompute, individualSpeakers=True)
            start = 0
            for i in range(len(micNumber)):
                nMic = micNumber[i]
                stop = nMic + start
                self.pMicArray.append(pMicArray[:, start:stop, :])
                self.pMic.append(pMic[:, start:stop])
                start = nMic + start

        elif type(obs_name) == list:
            nObs = len(obs_name)
            for i in range(nObs):
                tmpObs = obs_name[i]
                idxTmp = np.argwhere(np.array(self.observationName) == tmpObs)[0][0]
                isComputed = self.checkObservationExists(self.observationName[idxTmp])
                if isComputed is False:
                    micToCompute = np.concatenate((micToCompute, self.xMic[idxTmp]), axis=0)
                    self.computedObservations.append(self.observationName[idxTmp])
                    micNumber.append(len(self.xMic[idxTmp]))
                print('\n' + '=' * 40)
                # print('Solving pressure at: {}'.format(obs_name))
                pMic, pMicArray = self.bemObject.getMicPressure(micToCompute, individualSpeakers=True)
                start = 0
                for i in range(len(micNumber)):
                    nMic = micNumber[i]
                    stop = nMic + start
                    self.pMicArray.append(pMicArray[:, start:stop, :])
                    self.pMic.append(pMic[:, start:stop])
                    start = nMic + start

        elif type(obs_name) == str:
            tmpObs = obs_name
            idxTmp = np.argwhere(np.array(self.observationName) == tmpObs)[0][0]
            idxTmp = np.argwhere(np.array(self.observationName) == tmpObs)[0][0]
            isComputed = self.checkObservationExists(self.observationName[idxTmp])
            if isComputed is False:
                print('\n' + "=" * 40)
                # print('\n'+'Observation: {}'.format(self.observationName[idxTmp]))
                pMic, pMicArray = self.bemObject.getMicPressure(self.xMic[idxTmp], individualSpeakers=True)
                self.pMic.append(pMic)
                self.pMicArray.append(pMicArray)
                self.computedObservations.append(self.observationName[idxTmp])
        return None

    def deleteObservation(self, obsNameToDelete):
        if obsNameToDelete in self.observationName:
            index = self.observationName.index(obsNameToDelete)
            self.observationName.pop(index)
            self.observationType.pop(index)

            # xmic
            self.xMic.pop(index)
            self.labels.pop(index)

            # planar
            self.L.pop(index)
            self.W.pop(index)
            self.planarPlane.pop(index)

            # polar
            self.theta.pop(index)
            self.polarPlane.pop(index)
            self.DI = []

        if obsNameToDelete in self.computedObservations:
            index = self.computedObservations.index(obsNameToDelete)
            self.computedObservations.pop(index)
            self.pMic.pop(index)
            self.pMicArray.pop(index)
        return None

    def checkObservationExists(self, obsName, deleteObservation=False):
        """
        Check if observation exists in observationName or computedObservation.
        If already computed, nothing is done (microphone placement is not updated)
        If already exists but is not computed, the microphone placement is updated (observation is deleted then reset).
        If it does not exists, nothing is returned and the microphone placement is done as usual.
        :param obsName:
        :return:
        """
        if obsName in self.computedObservations:
            isComputed = True  # observation is already computed, no need to recompute it
        elif obsName in self.observationName:
            if deleteObservation is True:
                self.deleteObservation(obsName)
            isComputed = False # observation previously existed but wasn't computed, mic placement is updated
        else:
            isComputed = False
        return isComputed

    def plot(self, obs_name='all', radiatingSurface='all'):
        if obs_name == 'all':
            obs_name = self.computedObservations # return a list to compute everything in the next IF condition
        
        # reshape pMic or pMicArray in micData if not all radiating surface are plotted
        if radiatingSurface == 'all':
            micData = self.pMic
        elif radiatingSurface != 'all':
            micData = []
            for i in range(len(self.computedObservations)):
                ind_surface = []
                for j in range(len(radiatingSurface)):
                    ind_surface.append(np.argwhere(self.bemObject.radiatingSurface == radiatingSurface[j])[0][0])
                micData.append(np.sum(self.pMicArray[i][:, :, ind_surface], 2))

        toPlotter = []
        toPlotter_micIndex = []
        out = []
        if type(obs_name) == list:
            nObs = len(obs_name)
            for i in range(nObs):
                tmpObs = obs_name[i]
                pMicIdx = np.argwhere(np.array(self.computedObservations) == tmpObs)[0][0]
                idxTmp = np.argwhere(np.array(self.observationName) == tmpObs)[0][0]
                if self.observationType[idxTmp] == 'polar':
                    theta = np.array(self.theta[idxTmp], dtype=float)
                    cursor = plotTotalDirectivity(theta,
                                         self.bemObject.freq_array, micData[pMicIdx], title=self.observationName[idxTmp])
                    out.append(cursor)

                elif self.observationType[idxTmp] == 'spherical':
                    boxSphere = plot_spherical_radiation(self.xMic[idxTmp], micData[pMicIdx], self.bemObject.freq_array)
                    # out.append(boxSphere)

                elif self.observationType[idxTmp] == 'planar':
                    if bool(self.add2Plotter[idxTmp]) is False:
                        sfreq_spl = plotPressureField_SPL(self.L[idxTmp], self.W[idxTmp],
                                          micData[pMicIdx],
                                          self.bemObject.freq_array, 250,
                                          self.planarPlane[idxTmp], title=self.observationName[idxTmp])
                        sfreq_real = plotPressureField_REAL(self.L[idxTmp], self.W[idxTmp],
                                                      micData[pMicIdx],
                                                      self.bemObject.freq_array, 250,
                                                      self.planarPlane[idxTmp], title=self.observationName[idxTmp])
                        out.append(sfreq_spl)
                        out.append(sfreq_real)
                    elif bool(self.add2Plotter[idxTmp]) is True:
                        toPlotter.append(idxTmp)
                        toPlotter_micIndex.append(idxTmp)

                elif self.observationType[idxTmp] == 'FRF':
                    plotFRF(self.bemObject.freq_array, micData[pMicIdx],
                            self.labels[idxTmp])
                    out.append(None)

                elif self.observationType[idxTmp] == "BOX":
                    showPlot = plotBoundingBox_pyvista(self.bemObject,
                                                       self.Lx[idxTmp],
                                                       self.Ly[idxTmp],
                                                       self.Lz[idxTmp],
                                                       micData[pMicIdx], self.xMic[idxTmp], radiatingSurface)
                    
        elif type(obs_name) == str:
            tmpObs = obs_name
            pMicIdx = np.argwhere(np.array(self.computedObservations) == tmpObs)[0][0]
            idxTmp = np.argwhere(np.array(self.observationName) == tmpObs)[0][0]
            if self.observationType[idxTmp] == 'polar':
                theta = np.array(self.theta[idxTmp], dtype=float)
                cursor = plotTotalDirectivity(theta,
                                     self.bemObject.freq_array, micData[pMicIdx], title=self.observationName[idxTmp])
                out.append(cursor)

            elif self.observationType[idxTmp] == 'spherical':
                boxSphere = plot_spherical_radiation(self.xMic[idxTmp], micData[pMicIdx], self.bemObject.freq_array)
                # out.append(boxSphere)
            
            elif self.observationType[idxTmp] == 'planar':
                if self.add2Plotter[idxTmp] is False:
                    sfreq_spl = plotPressureField_SPL(self.L[idxTmp], self.W[idxTmp],
                                      micData[pMicIdx],
                                      self.bemObject.freq_array, 250,
                                      self.planarPlane[idxTmp], title=self.observationName[idxTmp])
                    sfreq_real = plotPressureField_REAL(self.L[idxTmp], self.W[idxTmp],
                                                  micData[pMicIdx],
                                                  self.bemObject.freq_array, 250,
                                                  self.planarPlane[idxTmp], title=self.observationName[idxTmp])
                    out.append(sfreq_spl)
                    out.append(sfreq_real)
                elif self.add2Plotter[idxTmp] is True:
                    toPlotter.append(idxTmp)
                    toPlotter_micIndex.append(pMicIdx)

            elif self.observationType[idxTmp] == 'FRF':
                plotFRF(self.bemObject.freq_array, micData[pMicIdx],
                        self.labels[idxTmp])
                out.append(None)

            elif self.observationType[idxTmp] == 'BOX':
                showPlot = plotBoundingBox_pyvista(self.bemObject,
                                                   self.Lx[idxTmp],
                                                   self.Ly[idxTmp],
                                                   self.Lz[idxTmp],
                                                   micData[pMicIdx], self.xMic[idxTmp], radiatingSurface)

        # plot pressure fields in PyVista
        if len(toPlotter) != 0:
            showPlot = plotPressureField_pyvista(self.bemObject, toPlotter,
                                                 toPlotter_micIndex, self.L,
                                                 self.W, micData, self.xMic, radiatingSurface)

        plt.show()
        return out


    ## POST PROCESSING FUNCTIONS
    def get_pMic(self, obsName, radiatingSurface='all', get_freq_array=False):
        """
        Return the pressure of the specified observation. pMic.shape = (freq_array, nMic)
        :param obsName: str
            observation name which will give the corresponding microphone pressure
        :return: out:
            pMic and frequency axis if get_freq_axis is True
        """

        # find corresponding observation
        pMicIdx = np.argwhere(np.array(self.computedObservations) == obsName)[0][0]

        # reshape pMic or pMicArray in micData if not all radiating surface are plotted
        if radiatingSurface == 'all':
            micData = self.pMic[pMicIdx]
            nMic = np.shape(self.pMic[pMicIdx])[1]
        elif radiatingSurface != 'all':
            ind_surface = []
            for j in range(len(radiatingSurface)):
                ind_surface.append(np.argwhere(self.bemObject.radiatingSurface == radiatingSurface[j])[0][0])
            micData = np.sum(self.pMicArray[pMicIdx][:, :, ind_surface], 2)
            nMic = np.shape(self.pMicArray[pMicIdx])[1]


        micData = micData.reshape([len(self.bemObject.freq_array), nMic])

        if get_freq_array is False:
            out = micData
        if get_freq_array is True:
            out = (micData, self.bemObject.freq_array)
        return out

    def get_directivityIndex(self, polarObsName):
        """
        Return directivity index of the given observation
        :param polarObsName: str
            observation on which to compute directivity index
        :return:
            show figure
        """

        # check if directivity index can be computed
        idxTmp = np.argwhere(np.array(self.observationName) == polarObsName)[0][0]
        pMicIdx = np.argwhere(np.array(self.computedObservations) == polarObsName)[0][0]  # polar obs must be computed
        if self.observationType[idxTmp] != 'polar':
            Exception('Observation must be polar.')


        # start and stop angles (index in array and associated value in radians)
        start_idx, start_val = gtb.findInArray(self.theta[idxTmp], 0)
        stop_idx, stop_val = gtb.findInArray(self.theta[idxTmp], np.pi)


        # compute Directivity factor
        Q = np.zeros([len(self.bemObject.freq_array)])
        pMicDir = self.pMic[pMicIdx]
        p_0 = pMicDir[:, start_idx]  # on-axis pressure

        for i in range(len(self.bemObject.freq_array)):
            numerator = (4*np.pi*np.abs((p_0/np.sqrt(2))**2))*(180/np.pi)
            sumPressure = 0
            for j in np.arange(start_idx, stop_idx+1, 1):
                sumPressure += np.abs((pMicDir[i, j]/np.sqrt(2))**2) * np.sin(self.theta[idxTmp][j])
            denominator = 2*np.pi*sumPressure
            Q[i] = numerator[i] / denominator

        fig, ax = plt.subplots()
        ax.semilogx(self.bemObject.freq_array, 10*np.log10(Q))
        ax.set(xlabel='Frequency [Hz]', ylabel='Directivity index [dB]')
        ax.grid(which='both', linestyle='dotted')
        plt.tight_layout()

        self.DI.append(10*np.log10(Q))
        return plt.show()


    def plot_system(self, engineering=False):
        mesh = pyvista.read(self.bemObject.meshPath)
        pl = pyvista.Plotter()
        pl.add_mesh(mesh, show_edges=True, cmap='summer', show_scalar_bar=False)
        light = pyvista.Light(light_type='headlight')
        # pl.add_floor('-z', color='grey', pad=5.0)

        pl.add_light(light)

        colour = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'yellow', 'cyan'] * 4

        for i in range(len(self.xMic)):
            point_cloud_tmp = pyvista.PolyData(self.xMic[i].astype(float))
            if self.observationType[i] == 'BOX':
                pl.add_mesh(point_cloud_tmp.outline(), color=colour[i],
                            render_points_as_spheres=True, label=self.observationName[i], point_size=9)
            else:
                pl.add_mesh(point_cloud_tmp, color=colour[i],
                            render_points_as_spheres=True, label=self.observationName[i], point_size=9)

        pl.add_axes(color='black')
        pl.add_legend(face='circle', bcolor=None)

        pl.background_color = 'white'

        _ = pl.show_grid(color='k')
        obj = pl.show()
        return obj


# =============================================================================
## MICROPHONE ARRAY FUNCTIONS
# =============================================================================
def createCircArray(theta, plane, radius=5, offset=[0,0,0]):
    """
    Create an array of shape (len(theta), 3) with the coordinates of a circular
    microphone array. To be used with pointSourceProblem.getMicPressure() or 
    bemppProblem.getMicPressure().

    Parameters
    ----------
    theta : numpy array
        Angles [radians] at which the observation points are placed.
    plane : str
        Plane on which the observation points are placed. ex: "xy", "zx", etc.
        First axis will define 0 degree.
    radius : float, optional
        Distance origin / microphones (radius of the circle) in meters. The 
        default is 5.
    offset : list of floats, optional
        Offset of the center of the circular array, in meters. 
        The default is [0,0,0].

    Raises
    ------
    ValueError
        Tell user if "plane" parameter is not understood - ex: "xa" will 
        raise this error.

    Returns
    -------
    numpy array
        mic_circle array of shape (len(theta), 3) with the x, y, z coordinates
        of a circular microphone array.

    """
    radius = radius
    xo = offset[0]
    yo = offset[1]
    zo = offset[2]

    if plane[0] == '-':
        sign = -1
        plane = plane[1:]
    else :
        sign = 1
        
    if plane == "xy":
        mic_circle = np.array([radius*np.cos(theta)+xo, 
                                radius*np.sin(theta)+yo, 
                                np.zeros(len(theta))+zo]) * np.array([[sign], [1], [1]])
    elif plane == "yx":
        mic_circle = np.array([radius*np.sin(theta)+xo, 
                               radius*np.cos(theta)+yo, 
                               np.zeros(len(theta))+zo]) * np.array([[1], [sign], [1]])
    elif plane == "xz":
        mic_circle = np.array([radius*np.cos(theta)+xo, 
                               np.zeros(len(theta))+yo, 
                               radius*np.sin(theta)+zo]) * np.array([[sign], [1], [1]])
    elif plane == "zx":
        mic_circle = np.array([radius*np.sin(theta)+xo, 
                               np.zeros(len(theta))+yo, 
                               radius*np.cos(theta)+zo]) * np.array([[1], [1], [sign]])
    elif plane == "yz":
        mic_circle = np.array([np.zeros(len(theta))+xo, 
                               radius*np.cos(theta)+yo,
                               radius*np.sin(theta)+zo]) * np.array([[1], [sign], [1]])
    elif plane == "zy":
        mic_circle = np.array([np.zeros(len(theta))+xo, 
                               radius*np.sin(theta)+yo,
                               radius*np.cos(theta)+zo]) * np.array([[1], [1], [sign]])
    else:
        raise ValueError("Plane not understood")
        
    return mic_circle.T


def createPlaneArray(length, width, micSpacing, plane, offset=[0, 0, 0], vert=False, mode=False):
    """
    Create a rectangular array of microphones on given plane. Place the corner on [x=0, y=0, z=0]

    Parameters
    ----------
    length : int or float
        rectangle length.
    width : int or float
        rectangle width.
    micSpacing : float
        distance from one microphone to another 
    plane : str
        plane 'xy', 'xz'. 'zy'.
    offset : list, optional
        Corner offset. The default is [0, 0, 0].

    Returns
    -------
    xmic : array
        Microphone location in a [nMic, 3] array. Compatible with 
        getMicPressure().
    L : array
        Array corresponding to the given length
    W : array
        Array corresponding to the given width

    """
    L = np.arange(0, length+micSpacing, micSpacing)
    W = np.arange(0, width+micSpacing, micSpacing)
    nMic = len(L)*len(W)
    xmic = np.zeros([nMic, 3])
    xOffset = np.ones([nMic, 3]) * offset
    # is there a better way to do it?
    if plane=='xy':
        dim1 = 0
        dim2 = 1
        dim3 = 2 
    elif plane=='yx':
        dim1 = 1
        dim2 = 0
        dim3 = 2
    elif plane=='xz':
        dim1 = 0
        dim2 = 2
        dim3 = 1
    elif plane=='zx':
        dim1 = 2
        dim2 = 0
        dim3 = 1
    elif plane=='yz':
        dim1 = 1
        dim2 = 2
        dim3 = 0
    elif plane=='zy':
        dim1 = 2
        dim2 = 1
        dim3 = 0
        
    i = 0
    for w in range(len(W)):
        for l in range(len(L)):
            xmic[i, dim1] = L[l]
            xmic[i, dim2] = W[w]
            xmic[i, dim3] = 0
            i += 1
    xmic += xOffset
    L += offset[dim1]
    W += offset[dim2]
    out = (xmic, L, W,)

    if vert is not False:
        xmic_n = filter_microphones(xmic, vert, mode=mode)
        L_n = L[np.isin(L, xmic_n[:, dim1])]
        W_n = W[np.isin(W, xmic_n[:, dim2])]
        out = (xmic_n, L_n, W_n)

    return out


def createSphereArray(nMic, sphereRadius=1.8, offset=[0, 0, 0]):
    """
    create a spherical microphone array
    :param Nmic: number of microphones in the array
    :return: xmic: cartesian coordinates of each microphone of the array
    """

    theta, phi = np.linspace(0, 2 * np.pi, int(np.sqrt(nMic))), np.linspace(0, np.pi, int(np.sqrt(nMic)))
    THETA, PHI = np.meshgrid(theta, phi)
    R = sphereRadius
    X = R * np.sin(PHI) * np.cos(THETA) + offset[0]
    Y = R * np.sin(PHI) * np.sin(THETA) + offset[1]
    Z = R * np.cos(PHI) + offset[2]
    
    mic = 0
    xmic = np.zeros([nMic, 3])
    for i in range(int(np.sqrt(nMic))):
        for j in range(int(np.sqrt(nMic))):
            xmic[mic, 0] = X[i, j]
            xmic[mic, 1] = Y[i, j]
            xmic[mic, 2] = Z[i, j]
            mic += 1

    return xmic

def create_boundingBox(Lx, Ly, Lz, step=1, offset=[0, 0, 0]):
    x_offset, y_offset, z_offset = offset
    x_range = np.arange(x_offset, x_offset + Lx+step, step)
    y_range = np.arange(y_offset, y_offset + Ly+step, step)
    z_range = np.arange(z_offset, z_offset + Lz+step, step)

    x_points, y_points, z_points = np.meshgrid(x_range, y_range, z_range)

    points = np.vstack([x_points.flatten(), y_points.flatten(), z_points.flatten()]).T

    return points


## PLOTTING FUNCTIONS
def plotFRF(freq_array, pmic, labels):
    nMic = np.shape(pmic)[1]
    fig, ax = plt.subplots()
    for i in range(nMic):
        ax.semilogx(freq_array, gtb.gain.SPL(pmic[:, i]), label=labels[i])
    ax.set(xlabel='Frequency [Hz]', ylabel='Sound Pressure Level [dB]')
    ax.grid(which='both')
    ax.legend(loc='best')
    return None

def plotTotalDirectivity(theta, freq, pMic, xscale='log', fmin=20, fmax=20e3,
                         dBmin=-20, dBmax=False, title="title"):
    """
    Plot the directivity from pMic obtained with getMicPressure for a circular microphone array

    Parameters
    ----------
    theta : narray
        observation angles.
    freq : narray
        frequency array.
    pMic : narray
        Pressur obtained with getMicPressure().
    xscale : str, optional
        scaling of the x axis (logarithmic, linear). The default is 'log'.
    fmin : float, optional
        min xlim. The default is 20.
    fmax : float, optional
        max xlim. The default is 20e3.
    dBmin : float, optional
        min ylim. The default is -20.
    dBmax : float, optional
        max ylim. The default is False.
    norm : bool, optional
        set the normalization. The default is True.

    Returns
    -------
    None.

    """
    print("THETA: ", theta)
    # theta = np.deg2rad(theta)
    # compute directivity
    maxAngle = np.array([np.max(abs(pMic), 1)])
    maxMatrix = np.repeat(maxAngle, len(theta), 0).T
    directivity = gtb.gain.dB(np.abs(pMic)/maxMatrix)
    dBmax = int(np.max(gtb.gain.SPL(pMic)))+3
    dBmin = dBmax-40
    
    fig = plt.figure(figsize=(11, 8))
    fig.suptitle(title)
    fig.subplots_adjust(top=0.88)
    ax1 = fig.add_subplot(221)
    gca1 = ax1.contourf(freq, np.rad2deg(theta), gtb.gain.SPL(pMic).T,
                             np.arange(dBmin, dBmax+3, 3), cmap='turbo')
    ax1.set_xscale('log')
    ax1.set(xlabel="Frequency [Hz]", ylabel="Angle [deg]", xlim=[fmin, fmax], 
            title='Raw directivity (SPL)')
    plt.colorbar(gca1)
    
    # Total directivity NORMALIZED
    ax2 = fig.add_subplot(223)
    gca2 = ax2.contourf(freq, np.rad2deg(theta), directivity.T, np.arange(-21, 3, 3), cmap='turbo')
    ax2.set_xscale('log')
    ax2.set(xlabel="Frequency [Hz]", ylabel="Angle [deg]", xlim=[fmin, fmax], 
            title='Absolute directivity')
    plt.colorbar(gca2)
    
    # pressure response
    obsAngle = 0
    ind_angle, value_angle = gtb.findInArray(theta, np.deg2rad(obsAngle))
    ax3 = fig.add_subplot(222)
    ax3.semilogx(freq, gtb.gain.SPL(pMic[:, ind_angle]), label=int(np.rad2deg(value_angle)))
    ax3.set(xlabel='Frequency [Hz]', ylabel='SPL [dB]', ylim=[dBmin, dBmax], 
            title='Pressure response')
    ax3.grid(which='both')
    ax3.legend(loc='best')
    
    # polar
    freqObs = 1000
    ind_freq, value_freq = gtb.findInArray(freq, freqObs)
    ax4 = fig.add_subplot(224, polar=True)
    ax4.plot(theta, gtb.gain.SPL(pMic[ind_freq, :]), label=int(value_freq))
    ax4.set(ylim=[dBmin, dBmax], title='Polar directivity')
    ax4.legend(loc='best')

    #defining the cursor
    cursor = Cursor(ax1, horizOn = True, vertOn=True, color='black', linewidth=1.2, 
                    useblit=True)
    
    def onclick(event):
        # global coord
        # coord.append((event.xdata, event.ydata))
        if event.button is MouseButton.LEFT:
            x = event.xdata
            y = event.ydata
             
            # update subplots -- FRF
            obsAngle = y
            ind_angle, value_angle = gtb.findInArray(theta, np.deg2rad(obsAngle))
            ax3.semilogx(freq, gtb.gain.SPL(pMic[:, ind_angle]), label=int(np.rad2deg(value_angle)))
            ax3.legend(loc='best')
            
            # update subplots -- polar
            freqObs = x
            ind_freq, value_freq = gtb.findInArray(freq, freqObs)
            ax4.plot(theta, gtb.gain.SPL(pMic[ind_freq, :]), label=int(value_freq))
            ax4.legend(loc='best')
            
            fig.canvas.draw() #redraw the figure
        
        if event.button is MouseButton.RIGHT:
            if len(ax3.lines) > 1:
                ax3.lines[-1].remove()
                ax4.lines[-1].remove()
                ax3.legend(loc='best')
                ax4.legend(loc='best')
                fig.canvas.draw()
        
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.show()
    return cursor

def plotPressureField_SPL(Length, Width, pMic,
                      freq_array, freq_obs, plane, normalization='SPL', 
                      dBLim=[50, 110], colorSteps=5, title="title"):
    """
    Plot 

    Parameters
    ----------
    Length : TYPE
        DESCRIPTION.
    Width : TYPE
        DESCRIPTION.
    xmic : TYPE
        DESCRIPTION.
    pMic : TYPE
        DESCRIPTION.
    freq_array : TYPE
        DESCRIPTION.
    freq_obs : TYPE
        DESCRIPTION.
    dBLim : TYPE, optional
        DESCRIPTION. The default is [50, 110].
    colorSteps : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # init data
    maxY = int(np.max(gtb.gain.SPL(pMic))) + 5
    minY = maxY - 75
    stepY = 5
    colorMap = 'turbo'

    init_frequency = 150
    fmin = freq_array[0] 
    fmax = freq_array[-1]
    dim1 = plane[0]
    dim2 = plane[1]

    # processed output values
    f_index, f_value = gtb.findInArray(freq_array, freq_obs)
    pMicRes = np.reshape(pMic[f_index, :], [len(Width), len(Length)])
    pMicSPL = gtb.gain.SPL(pMic)
    pMicReal = np.real(pMic)
    results = pMicSPL

    # initial plot
    fig, ax = plt.subplots()
    gca = ax.contourf(Length, Width, gtb.gain.SPL(pMicRes), np.arange(minY, maxY, stepY), cmap=colorMap)
    plt.subplots_adjust(bottom=0.25)
    plt.colorbar(gca, label='SPL [dB]')
    ax.axis('equal')
    ax.set(xlabel= dim1 + ' [m]', ylabel= dim2 + ' [m]', title=title + ' - f = '+ str(int(f_value)) + ' Hz')

    def update(expr):
        ax.clear()
        f_index, f_value = gtb.findInArray(freq_array, int(expr))
        pMicRes = np.reshape(results[f_index, :], [len(Width), len(Length)])
        
        ax.contourf(Length, Width, pMicRes,
                          np.arange(minY, maxY, stepY), cmap=colorMap)
        ax.set_title("f = {} Hz".format(int(f_value)))
        ax.set(xlabel=dim1 + ' [m]', ylabel=dim2 + ' [m]', title=title + ' - f = '+ str(int(f_value)) + ' Hz')
        plt.draw()



    # sfreq
    graphBox = fig.add_axes([0.2, 0.05, 0.3, 0.075])
    sfreq = TextBox(graphBox, "f (Hz): ")
    sfreq.on_submit(update)
    sfreq.set_val(str(freq_obs))
    return sfreq


def plotPressureField_REAL(Length, Width, pMic,
                         freq_array, freq_obs, plane, normalization='SPL',
                         dBLim=[50, 110], colorSteps=5, title="title"):
    """
    Plot

    Parameters
    ----------
    Length : TYPE
        DESCRIPTION.
    Width : TYPE
        DESCRIPTION.
    xmic : TYPE
        DESCRIPTION.
    pMic : TYPE
        DESCRIPTION.
    freq_array : TYPE
        DESCRIPTION.
    freq_obs : TYPE
        DESCRIPTION.
    dBLim : TYPE, optional
        DESCRIPTION. The default is [50, 110].
    colorSteps : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # init data
    maxY = 1+0.1
    minY = -1
    stepY = 0.1
    colorMap = 'seismic'

    init_frequency = 150
    fmin = freq_array[0]
    fmax = freq_array[-1]
    dim1 = plane[0]
    dim2 = plane[1]

    # processed output values
    f_index, f_value = gtb.findInArray(freq_array, freq_obs)
    pMicReal = np.real(pMic)
    results = pMicReal
    pMicRes = np.reshape(results[f_index, :], [len(Width), len(Length)])

    # initial plot
    fig, ax = plt.subplots()
    gca = ax.contourf(Length, Width, pMicRes, np.arange(minY, maxY, stepY), cmap=colorMap)
    plt.subplots_adjust(bottom=0.25)
    plt.colorbar(gca, label='Normalized real pressure')
    ax.axis('equal')
    ax.set(xlabel=dim1 + ' [m]', ylabel=dim2 + ' [m]', title=title + ' - f = ' + str(int(f_value)) + ' Hz')

    def update(expr):
        ax.clear()
        f_index, f_value = gtb.findInArray(freq_array, int(expr))
        pMicRes = np.reshape(results[f_index, :], [len(Width), len(Length)])

        ax.contourf(Length, Width, pMicRes,
                    np.arange(minY, maxY, stepY), cmap=colorMap)
        ax.set_title("f = {} Hz".format(int(f_value)))
        ax.set(xlabel=dim1 + ' [m]', ylabel=dim2 + ' [m]', title=title + ' - f = ' + str(int(f_value)) + ' Hz')
        plt.draw()

    # sfreq
    graphBox = fig.add_axes([0.2, 0.05, 0.3, 0.075])
    sfreq = TextBox(graphBox, "f (Hz): ")
    sfreq.on_submit(update)
    sfreq.set_val(str(freq_obs))
    return sfreq

def plotPressureField_pyvista(bemOBJ, toPlotter, toPlotterMicIdx, L, W, pMicData, xMic, radiatingSurface):
    """
    Plot pressure field using PyVista backend.
    :param bemOBJ: object,
        bem object used in the current observation.
    :param toPlotter: list,
        list of observations to add in PyVista Plotter.
    :param toPlotterMicIdx: list,
        list of computed observations to add in PyVista Plotter.
    :param L: list,
        list of length of each planar surface defined using "planar" observations.
    :param W: list,
        list of width of each planar surface defined using "planar" observations.
    :param pMicData: ndarray,
        pressure field data.
    :param xMic: ndarray,
        microphones position.
    :param radiatingSurface: str, list,
        Index of surfaces currently radiating.
    :return:
    """
    freq_array = bemOBJ.freq_array
    meshPath = bemOBJ.meshPath
    meshPlot = []

    # init plotter
    pl = pyvista.Plotter(lighting='three lights')
    pl.background_color = 'white'
    sizeFactor = bemOBJ.sizeFactor
    vertices = bemOBJ.vertices

    for i, idxData in enumerate(toPlotterMicIdx):
        meshPlot.append(pyvista.StructuredGrid())
        # initialisation
        data = pMicData[idxData].T
        meshPlot[i].points = xMic[idxData]
        meshPlot[i].points = xMic[idxData]
        meshPlot[i].dimensions = [len(L[idxData]), len(W[idxData]), 1]


    engine = pyvista_spkPlotter(pl, bemOBJ, toPlotter, toPlotterMicIdx, L, W,
                                pMicData, xMic, radiatingSurface)
    # chart = pyvista.Chart2D()

    def plotPressurePoint(point):
        position = np.zeros([len(toPlotter), 3])
        for i in range(len(toPlotter)):
            _, position[i, :] = gtb.geometry.findClosestPoint(xMic[toPlotterMicIdx[i]], point[0], point[1], point[2])
        # print(position[i, :])
        idxp = gtb.geometry.findClosestPoint(position, point[0], point[1], point[2])[0]
        idxm, p = gtb.geometry.findClosestPoint(xMic[toPlotterMicIdx[idxp]], point[0], point[1], point[2])

        pMic_tmp = pMicData[toPlotterMicIdx[idxp]].T
        pMicToPlot = pMic_tmp[idxm, :]
        max = np.max(gtb.gain.SPL(pMicToPlot))

        fig, ax = plt.subplots()
        ax.semilogx(bemOBJ.freq_array, gtb.gain.SPL(pMicToPlot), label='{}'.format(p))
        ax.grid(which='both')
        ax.set(xlabel='Frequency [Hz]', ylabel='SPL [dB]', ylim=[max-40, max+6])
        ax.legend(loc='best')
        plt.show()

        # chart = pyvista.Chart2D()
        # chart.line(np.log10(bemOBJ.freq_array), gtb.gain.SPL(pMicToPlot), label='{}'.format(point))
        # chart.x_label='Frequency [kHz]'
        # chart.y_label='SPL [dB]'
        # chart.show()
        return None


    pl.add_slider_widget(
        callback=lambda value: engine('cfreq', value),
        rng=[np.log10(freq_array[0]), np.log10(freq_array[-1])],
        value=np.log10(250),
        title="Frequency [log]",
        style="modern"
    )

    pl.add_checkbox_button_widget(
        callback=lambda value: engine('real', value),
        value=False,
        position=(5, 100),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('showMesh', value),
        value=False,
        position=(5, 0),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('contourP', value),
        value=True,
        position=(5, 50),
    )

    pl.enable_surface_point_picking(
        callback=plotPressurePoint,
        show_message=False,
    )

    obj = pl.show()
    return obj


def plotBoundingBox_pyvista(bemOBJ, nx, ny, nz, pMicData, xMic, radiatingSurface):
    """
    Plot a bounding box observation in a PyVista Plotter
    :param bemOBJ: object,
        bem object used in the current observation.
    :param pMicData: ndarray,
        pressure field data.
    :param xMic: ndarray,
        microphones position.
    :param radiatingSurface: str, list,
        Index of surfaces currently radiating.
    :return:
    """
    freq_array = bemOBJ.freq_array
    meshPath = bemOBJ.meshPath
    meshPlot = []

    # init plotter
    pl = pyvista.Plotter(lighting='three lights')
    pl.background_color = 'white'
    sizeFactor = bemOBJ.sizeFactor
    vertices = bemOBJ.vertices

    engine = pyvista_boundingBoxPlotter(pl, bemOBJ, nx, ny, nz, pMicData, xMic, radiatingSurface)

    pl.add_slider_widget(
        callback=lambda value:engine('cfreq', value),
        rng=[np.log10(freq_array[0]), np.log10(freq_array[-1])],
        value=np.log10(250),
        title="Frequency, log()",
        style="modern"
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('real', value),
        value=False,
        position=(5, 100),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('showMesh', value),
        value=False,
        position=(5, 0),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('contourP', value),
        value=True,
        position=(5, 50),
    )

    pl.add_plane_widget(
        callback=lambda normal, origin: engine('normOrigin', (normal, origin)),
        normal=(1, 0, 0),
        origin=(0, 0, 0),
        normal_rotation=True,
    )

    obj = pl.show()
    return obj

def plotSphericalRadiation(xMic, pMic, freq_array):
    """
    Balloon plot of pressure on spherical microphone array
    :param xMic:
    :param pMic:
    :param freq_array:
    :return:
    """

    Nfft = len(freq_array)
    nMic = len(xMic)
    nMicSq = int(np.sqrt(nMic))
    X = np.zeros([nMicSq, nMicSq, Nfft])
    Y = np.zeros([nMicSq, nMicSq, Nfft])
    Z = np.zeros([nMicSq, nMicSq, Nfft])


    # idxFreq, valueFreq = gtb.findInArray(freq_array, 100)
    spl = gtb.gain.dBSPL(pMic)
    SPL = np.zeros([nMicSq, nMicSq, Nfft])
    mic = 0
    for i in range(nMicSq):
        for j in range(nMicSq):
            for k in range(Nfft):
                X[i, j, k] = xMic[mic, 0] / np.max(np.abs(xMic)) * spl[k, mic]
                Y[i, j, k] = xMic[mic, 1] / np.max(np.abs(xMic)) * spl[k, mic]
                Z[i, j, k] = xMic[mic, 2] / np.max(np.abs(xMic)) * spl[k, mic]
                SPL[i, j, k] = spl[k, mic]
            mic += 1

    # SPL value for the balloon plot
    idxFreq, valueFreq = gtb.findInArray(freq_array, 100)
    cmap = plt.get_cmap('turbo')
    norm = colors.Normalize(vmin=SPL[:, :, idxFreq].min(), vmax=SPL[:, :, idxFreq].max())

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    plot = ax.plot_surface(
        X[:, :, idxFreq], Y[:, :, idxFreq], Z[:, :, idxFreq], rstride=1, cstride=1,
        facecolors=cmap(norm(SPL[:, :, idxFreq])),
        linewidth=0, antialiased=False, alpha=1, shade=False)
    ax.axis('equal')
    ax.set(xlabel='x [dB]', ylabel="y [dB]", zlabel="z [dB]", title="Acoustic projection at F={}Hz".format(int(valueFreq)))

    # data update
    def update(expr):
        ax.clear()
        idxFreq, valueFreq = gtb.findInArray(freq_array, int(expr))
        cmap = plt.get_cmap('turbo')
        norm = colors.Normalize(vmin=SPL[:, :, idxFreq].min(), vmax=SPL[:, :, idxFreq].max())

        plot = ax.plot_surface(
            X[:, :, idxFreq], Y[:, :, idxFreq], Z[:, :, idxFreq], rstride=1, cstride=1,
            facecolors=cmap(norm(SPL[:, :, idxFreq])),
            linewidth=0, antialiased=False, alpha=1, shade=False)
        ax.axis('equal')
        ax.set(xlabel='x [dB]', ylabel="y [dB]", zlabel="z [dB]",
               title="Acoustic projection at F={}Hz".format(int(valueFreq)))
        plt.draw()


    # sfreq
    graphBox = fig.add_axes([0.2, 0.05, 0.3, 0.075])
    sfreq = TextBox(graphBox, "f (Hz): ")
    sfreq.on_submit(update)
    sfreq.set_val(str(100))

    return sfreq


def plot_spherical_radiation(xMic, pMic, freq_array):
    # init plotter
    pl = pyvista.Plotter(lighting='three lights')
    pl.background_color = 'white'

    engine = pyvista_sphericalPlotter(pl, pMic, xMic, freq_array)

    pl.add_slider_widget(
        callback=lambda value:engine('cfreq', value),
        rng=[np.log10(freq_array[0]), np.log10(freq_array[-1])],
        value=np.log10(250),
        title="Frequency, log()",
        style="modern"
    )
    obj = pl.show()
    return obj


## tools
def sumPressureArray(bemObj, radiatingSurface):
    # print(radiatingSurface)
    p_mesh = bemObj.p_mesh
    radSurf_system = bemObj.radiatingSurface
    if radiatingSurface == 'all':
        radiatingSurface = radSurf_system

    pressureCoeff = np.zeros([len(bemObj.freq_array), bemObj.spaceP.grid_dof_count], dtype='complex')
    for f in range(len(bemObj.freq_array)):
        for j in range(len(radiatingSurface)):
            ind_surface = np.argwhere(radSurf_system == radiatingSurface[j])[0][0]
            pressureCoeff[f, :] += p_mesh[f][ind_surface].coefficients
    return pressureCoeff


def filter_microphones(microphones, boundary, mode='inside'):
    """
    Filter microphones based on their location relative to the boundary.

    Parameters:
    - microphones: NumPy array of shape (N, 3), representing microphone coordinates.
    - boundary: NumPy array of shape (M, 3), representing boundary coordinates.
    - mode: 'inside' or 'outside', specifying whether to keep microphones inside or outside the boundary.

    Returns:
    - filtered_microphones: NumPy array of shape (K, 3), where K <= N, containing the filtered microphone coordinates.
    """

    if mode not in ('inside', 'outside'):
        raise ValueError("Mode must be 'inside' or 'outside'.")

    if mode == 'inside':
        # Keep microphones inside the boundary
        condition = np.all((microphones[:, np.newaxis] >= boundary.min(axis=0)) &
                           (microphones[:, np.newaxis] <= boundary.max(axis=0)), axis=2)
    else:
        # Keep microphones outside the boundary
        condition = np.any((microphones[:, np.newaxis] < boundary.min(axis=0)) |
                           (microphones[:, np.newaxis] > boundary.max(axis=0)), axis=2)

    filtered_microphones = microphones[condition.all(axis=1)]
    return filtered_microphones



## ================================================
# %% PyVista plotter classes (for 3D visualization)
class pyvista_spkPlotter:
    def __init__(self, mesh, bemOBJ, toPlotter, toPlotterMicIdx, L, W, pMicData, xMic, radiatingSurface):
        self.output = mesh
        self.bemOBJ = bemOBJ
        self.freq_array = bemOBJ.freq_array
        self.vertices = bemOBJ.vertices
        self.sizeFactor = bemOBJ.sizeFactor
        self.toPlotter = toPlotter
        self.toPlotterMicIdx = toPlotterMicIdx
        self.L = L
        self.W = W
        self.pMicData = pMicData
        self.xMic = xMic
        self.radiatingSurface = radiatingSurface

        self.speaker = pyvista.read(bemOBJ.meshPath)
        meshPlot = []
        # pmicdata
        for i, idxData in enumerate(toPlotterMicIdx):
            meshPlot.append(pyvista.StructuredGrid())
            # initialisation
            data = pMicData[idxData].T
            meshPlot[i].points = xMic[idxData]
            # pMicRes = np.reshape(data[:, 0], [len(L[idxData]), len(W[idxData])])
            meshPlot[i].points = xMic[idxData]
            meshPlot[i].dimensions = [len(L[idxData]), len(W[idxData]), 1]
        self.meshPlot = meshPlot


        # default parameters
        self.kwargs = {
            'real': False,
            'cfreq': '250',
            'showMesh': False,
            'contourP': True
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        bemOBJ = self.bemOBJ
        real = self.kwargs['real']
        cfreq = self.kwargs['cfreq']
        showMesh = self.kwargs['showMesh']
        contourP = self.kwargs['contourP']


        # min / max bounds - title
        maxPressure = []
        minPressure = []
        if real is True:
            for i in self.toPlotterMicIdx:
                maxPressure.append(np.max(np.real(self.pMicData[i])))
                minPressure.append(np.min(np.real(self.pMicData[i])))
            minPressure = np.min(minPressure)
            maxPressure = np.max(maxPressure)
            pressureCoeff_system = sumPressureArray(self.bemOBJ, self.radiatingSurface)
            pressureCoeff_system = gtb.normMinMax(np.real(pressureCoeff_system), minPressure, maxPressure)
            title_arg = 'Real(Pressure), Pa - '
            cmap_d = 'seismic'
        else:
            for i in self.toPlotterMicIdx:
                maxPressure.append(np.max(gtb.gain.SPL(self.pMicData[i])))

            maxPressure = np.max(maxPressure) + 5
            minPressure = maxPressure - 75

            # to normalize pressureCoeff_system
            # (pressure on baffle will usually be super high compared to radiated field)
            maxP_abs = 10**(maxPressure/20) * 2e-5

            pressureCoeff_system = sumPressureArray(self.bemOBJ, self.radiatingSurface)
            pressureCoeff_system = gtb.normMinMax(np.abs(pressureCoeff_system), 0, maxP_abs)
            pressureCoeff_system = gtb.gain.SPL(pressureCoeff_system)
            title_arg = 'SPL, dB - '
            cmap_d = 'turbo'

        if contourP is True:
            ncolor=21
        else:
            ncolor=256

        freq_idx, frequency = gtb.findInArray(self.freq_array, 10 ** (cfreq))
        sargs = dict(
            title_font_size=21,
            label_font_size=18,
            shadow=False,
            n_labels=5,
            italic=False,
            font_family="arial",
            color='k',
            title=title_arg + str(round(frequency, 1)) + " Hz",
            # position_x=0.35,
            # position_y=0
        )

        show_sc_b = []
        for i, idxData in enumerate(self.toPlotterMicIdx):
            if i == 0:
                show_sc_b.append(True)
            else:
                show_sc_b.append(False)

            data = self.pMicData[idxData].T
            pMicRes = np.reshape(data[:, int(freq_idx)],
                                 [len(self.L[idxData]), len(self.W[idxData])])

            # def real or SPL pressure plot
            if real is True:
                pMicToPlot = gtb.normMinMax(np.real(pMicRes), minPressure, maxPressure)
            else:
                pMicToPlot = gtb.gain.SPL(pMicRes)

            # system mesh
            self.output.remove_actor('spk_mesh')
            try:
                self.output.add_mesh(self.speaker, show_edges=showMesh,
                            scalars=pressureCoeff_system[freq_idx, :self.bemOBJ.spaceP.grid_dof_count // self.sizeFactor],
                            cmap=cmap_d,
                            n_colors=ncolor,
                            show_scalar_bar=False, clim=[minPressure, maxPressure], name='spk_mesh')
            except:
                pressure_system = pressureCoeff_system[freq_idx, :self.bemOBJ.spaceP.grid_dof_count // self.sizeFactor]
                toAdd = self.vertices - self.bemOBJ.spaceP.grid_dof_count
                # print("ZERO PADDING:::: ", toAdd)
                scalars_to_plot = np.concatenate((pressure_system, np.zeros(toAdd)))
                self.output.add_mesh(self.speaker, show_edges=showMesh,
                                     scalars=scalars_to_plot,
                                     cmap=cmap_d,
                                     n_colors=ncolor,
                                     show_scalar_bar=False, clim=[minPressure, maxPressure],
                                     name='spk_mesh')

            # observation mesh
            self.output.remove_actor('SPL_{}'.format(i))
            self.output.add_mesh(self.meshPlot[i], show_edges=showMesh, scalars=pMicToPlot, cmap=cmap_d,
                        scalar_bar_args=sargs, clim=[minPressure, maxPressure], n_colors=ncolor,
                        name='SPL_{}'.format(i), show_scalar_bar=show_sc_b[i])
        _ = self.output.show_grid(color='k')
        self.output.add_axes(color='k')
        return


class pyvista_boundingBoxPlotter:
    def __init__(self, mesh, bemOBJ, nx, ny, nz, pMicData, xMic, radiatingSurface):
        self.output = mesh
        self.bemOBJ = bemOBJ
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.freq_array = bemOBJ.freq_array
        self.vertices = bemOBJ.vertices
        self.sizeFactor = bemOBJ.sizeFactor
        self.pMicData = pMicData
        # self.pMicDataCopy = copy(self.pMicData)
        self.xMic = xMic
        self.radiatingSurface = radiatingSurface

        self.speaker = pyvista.read(bemOBJ.meshPath) # will display the speaker's mesh

        # pmicdata
        meshPlot = pyvista.StructuredGrid()
        meshPlot.points = xMic
        meshPlot.dimensions = [nx, ny, nz]
        self.meshPlot = meshPlot # will display bounding box mesh (unstructured grid?)

        # default parameters
        self.kwargs = {
            'real': False,
            'cfreq': '250',
            'showMesh': False,
            'contourP': True,
            'normal': (1, 0, 0),
            'origin': (0, 0, 0)
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        bemOBJ = self.bemOBJ
        real = self.kwargs['real']
        cfreq = self.kwargs['cfreq']
        showMesh = self.kwargs['showMesh']
        contourP = self.kwargs['contourP']
        normOrigin = self.kwargs['normOrigin']
         #origin = self.kwargs['origin']

        normal = normOrigin[0]
        origin = normOrigin[1]
        # print(normal)
        # print(origin)
        # min / max bounds - title
        if real is True:
            minPressure = np.min(np.real(self.pMicData))
            maxPressure = np.max(np.real(self.pMicData))
            pressureCoeff_system = sumPressureArray(self.bemOBJ, self.radiatingSurface)
            pressureCoeff_system = gtb.normMinMax(np.real(pressureCoeff_system), minPressure, maxPressure)
            title_arg = 'Real(Pressure), Pa - '
            cmap_d = 'seismic'
        else:
            maxPressure = np.max(gtb.gain.SPL(self.pMicData)) + 5
            minPressure = maxPressure - 75

            # to normalize pressureCoeff_system
            # (pressure on baffle will usually be super high compared to radiated field)
            maxP_abs = 10 ** (maxPressure / 20) * 2e-5

            pressureCoeff_system = sumPressureArray(self.bemOBJ, self.radiatingSurface)
            pressureCoeff_system = gtb.normMinMax(np.abs(pressureCoeff_system), 0, maxP_abs)
            pressureCoeff_system = gtb.gain.SPL(pressureCoeff_system)
            title_arg = 'SPL, dB - '
            cmap_d = 'turbo'

        if contourP is True:
            ncolor = 21
        else:
            ncolor = 256

        freq_idx, frequency = gtb.findInArray(self.freq_array, 10 ** (cfreq))
        sargs = dict(
            title_font_size=21,
            label_font_size=18,
            shadow=False,
            n_labels=5,
            italic=False,
            font_family="arial",
            color='k',
            title=title_arg + str(round(frequency, 1)) + " Hz",
        )


        data = self.pMicData.T
        pMicRes = data[:, int(freq_idx)]

        # def real or SPL pressure plot
        if real is True:
            pMicToPlot = gtb.normMinMax(np.real(pMicRes), minPressure, maxPressure)
        else:
            pMicToPlot = gtb.gain.SPL(pMicRes)

        # system mesh
        _ = self.output.remove_actor('spk_mesh')
        self.output.add_mesh(self.speaker, show_edges=showMesh,
                             scalars=pressureCoeff_system[freq_idx, :self.vertices // self.sizeFactor],
                             cmap=cmap_d, n_colors=ncolor, show_scalar_bar=False,
                             clim=[minPressure, maxPressure], name='spk_mesh')

        # observation mesh
        meshPlot = self.meshPlot
        meshPlot.point_data['scalars'] = pMicToPlot
        slc = meshPlot.slice(normal=normal, origin=origin)
        # self.output.add_mesh(meshPlot.outline(), color='k')
        self.output.add_mesh(slc, show_edges=showMesh, show_scalar_bar=True,
                                   cmap=cmap_d,  clim=[minPressure, maxPressure], n_colors=ncolor,
                                   name='SPL_slice', scalar_bar_args=sargs)
        _ = self.output.show_grid(color='k')
        self.output.add_axes(color='k')
        return


class pyvista_sphericalPlotter:
    def __init__(self, mesh, pMic, xMic, frequency):
        self.output = mesh
        self.pMic = pMic
        self.xMic = xMic
        self.frequency = frequency

        ELEVATION = np.arccos(xMic[:, 2] / np.sqrt(xMic[:, 0] ** 2 + xMic[:, 1] ** 2 + xMic[:, 2] ** 2))
        AZIMUTH = np.zeros(len(xMic))
        self.R = np.zeros([len(frequency), len(xMic)])
        for i in range(len(xMic)):
            if xMic[i, 1] > 0:
                AZIMUTH[i] = np.arccos(xMic[i, 0] / np.sqrt(xMic[i, 0] ** 2 + xMic[i, 1] ** 2))
            elif xMic[i, 1] < 0:
                AZIMUTH[i] = - np.arccos(xMic[i, 0] / np.sqrt(xMic[i, 0] ** 2 + xMic[i, 1] ** 2))
            # self.R = gtb.gain.SPL(pMic[])

        self.elevation = ELEVATION
        self.azimuth = AZIMUTH

        self.kwargs = {
            'cfreq': '250',
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        cfreq = self.kwargs['cfreq']
        # print(cfreq)
        R = self.R
        elevation = self.elevation
        azimuth = self.azimuth

        freq_idx, frequency = gtb.findInArray(self.frequency, 10 ** (cfreq))
        # print(self.pMic.shape)

        sargs = dict(
            title_font_size=21,
            label_font_size=18,
            shadow=False,
            n_labels=5,
            italic=False,
            font_family="arial",
            color='k',
            title="SPL, dB - " + str(round(frequency, 1)) + " Hz",
        )

        # X = gtb.gain.SPL(np.sin(elevation) * np.cos(azimuth) * self.pMic[freq_idx, :]) #* np.sin(elevation) * np.cos(azimuth)
        # Y = gtb.gain.SPL(np.sin(elevation) * np.sin(azimuth) * self.pMic[freq_idx, :]) #* np.sin(elevation) * np.sin(azimuth)
        # Z = gtb.gain.SPL(np.cos(elevation) * self.pMic[freq_idx, :]) #* np.cos(elevation)
        # XVEC = np.zeros([len(X), 3])
        # XVEC[:, 0] = X
        # XVEC[:, 1] = Y
        # XVEC[:, 2] = Z
        cloud = pyvista.PolyData(self.xMic)
        volume = cloud.delaunay_3d(alpha=2.)
        shell = volume.extract_geometry()
        shell.points /= np.sqrt(self.xMic[0, 0]**2 + self.xMic[0, 1]**2 + self.xMic[0, 2]**2)
        # print(np.sqrt(self.xMic[0, 0]**2 + self.xMic[0, 1]**2 + self.xMic[0, 2]**2))
        shell.points *= np.tile(gtb.gain.SPL(self.pMic[freq_idx, :]), (3, 1)).T
        # points = np.copy(shell.points)

        _ = self.output.remove_actor('balloon')
        self.output.add_mesh(shell, name='balloon', cmap="turbo",
                             scalars=gtb.gain.SPL(self.pMic[freq_idx, :]),
                             clim=[np.max(gtb.gain.SPL(self.pMic))-30, np.max(gtb.gain.SPL(self.pMic)) + 6],
                             n_colors=20, scalar_bar_args=sargs)
        _ = self.output.show_grid(color='k')
        self.output.add_axes(color='k')
        return None

