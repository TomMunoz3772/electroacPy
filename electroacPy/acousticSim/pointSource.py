#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:32:22 2025

@author: tom
"""

import bempp.api
from bempp.api.operators.boundary import helmholtz, sparse
from bempp.api.operators.potential import helmholtz as helmholtz_potential
from bempp.api.assembly.discrete_boundary_operator import DiagonalOperator
from scipy.sparse.linalg import gmres as scipy_gmres
from bempp.api.linalg import gmres
import numpy as np
from tqdm import tqdm
import warnings
import electroacPy.general as gtb
from .ACSHelpers_ import (getSurfaceAdmittance, mirror_mesh, 
                         greenMonopole, incidenceCoeff, buildGridFunction,
                         element2vertice_pressure)

warnings.filterwarnings("ignore", message="splu requires CSC matrix format")
warnings.filterwarnings("ignore", message="splu converted its input to CSC format")

try:
    from pyopencl import CompilerWarning
    warnings.filterwarnings("ignore", category=CompilerWarning)
except:
    None


class pointSource:
    def __init__(self, sourcePosition, volumeVelocity, frequency, c=343, 
                 rho=1.22, boundary_conditions=None, **kwargs):
        """
        Simulate point sources in free space.

        Parameters
        ----------
        sourcePosition : list of list
            Source(s) coordinates [[x1, y1, z1], [x2, y2, z2], ..., [xN, yN, zN]].
        volumeVelocity : list of numpy array
            Sources volume velocity.
        frequency : numpy array
            Range of study.
        c : float, optional
            Speed of sound in propagation medium. The default is 343 m/s (air).
        rho : float, optional
            Density of propagation medium. The default is 1.22 kg/m^3.
        boundary_conditions : boundaryConditions object, optional
            Set of conditions on boundaries / infinite boundaries. The default is None.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        self.xSource = np.array(sourcePosition)
        self.xSystem = None
        self.QSource = np.array(volumeVelocity)
        self.QSystem = None
        self.frequency = frequency
        self.k = self.frequency*2*np.pi / c
        self.c = c
        self.rho = rho
        self.boundary_conditions = boundary_conditions
        self.sizeFactor = 1
        
        self.Nsource = len(self.xSource)
        if self.boundary_conditions is not None:
            self.mirror_system()
        else:
            self.xSystem = self.xSource
            
        self.isComputed = False

        # compatibility layer for integration with already-existing BEM eval 
        # and plots        
        if "radiatingElement" in kwargs:
            self.radiatingElement = kwargs["radiatingElement"]
        self.meshPath = None


    def mirror_system(self):
        self.xSystem = np.copy(self.xSource)
        for param in self.boundary_conditions.parameters:
            if param in ["x", "X"]:
                xmirror = np.copy(self.xSystem)
                xmirror[:, 0] = 2*self.boundary_conditions.parameters[param]["offset"] - xmirror[:, 0]
                self.xSystem = np.concatenate([self.xSystem, xmirror])
                self.sizeFactor *= 2
            if param in ["y", "Y"]:
                ymirror = np.copy(self.xSystem)
                ymirror[:, 1] = 2*self.boundary_conditions.parameters[param]["offset"] - ymirror[:, 1]
                self.xSystem = np.concatenate([self.xSystem, ymirror])
                self.sizeFactor *= 2
            if param in ["z", "Z"]:
                zmirror = np.copy(self.xSystem)
                zmirror[:, 2] = 2*self.boundary_conditions.parameters[param]["offset"] - zmirror[:, 2]
                self.xSystem = np.concatenate([self.xSystem, zmirror])
                self.sizeFactor *= 2
        self.QSystem = np.tile(self.QSource, (self.sizeFactor, 1))

    def getMicPressure(self, xMic, individualSpeakers=True):
        """
        Compute pressure at evaluation points (xMic).

        Parameters
        ----------
        xMic : numpy array
            Evaluation points in 3D space: shape (Nmic, 3).

        Returns
        -------
        pmic : numpy array
            Evaluated pressure (complex). Shape (Nfft, Nmic, Nsource)

        """
        pmic = np.zeros([len(self.frequency), len(xMic), len(self.xSystem)],
                        dtype=complex) 
        for m, mic in enumerate(xMic):
            for s, source in enumerate(self.xSystem):
                R = ((source[0]-mic[0])**2 + 
                    (source[1]-mic[1])**2 + 
                    (source[2]-mic[2])**2)**0.5
                pmic[:, m, s] = (1j*self.k*self.rho*self.c*self.QSystem[s, :] * 
                                 np.exp(-1j*self.k*R) / (4*np.pi*R))
        
        # need to sum all corresponding image/source group
        if self.sizeFactor != 1:
            Nsources = len(self.xSource)
            Nimages = self.sizeFactor * Nsources - Nsources 
            for s in range(Nsources):
                for i in range(Nimages//2):
                    pmic[:, :, s] += pmic[:, :, s + Nsources*i]
                            
        pmic = pmic[:, :, :Nsources]
        
        if individualSpeakers is True:
            out = (np.sum(pmic, 2), pmic)
        elif individualSpeakers is False:
            out = np.sum(pmic, 2)
        return out


class pointSourceBEM:
    def __init__(self, meshPath, xSource, volumeVelocity, frequency, 
             domain="interior", c_0=343, rho_0=1.22, **kwargs):
        """
        Simulate point source with additional boundary-element mesh for diffraction. 
        Mostly useful for room acoustics.

        Parameters
        ----------
        meshPath : str
            Path to simulation mesh.
        xSource : list of list
            Source(s) coordinates.
        velocity : list of numpy array
            Source volume velocity.
        frequency : numpy array
            Range of study.
        domain : TYPE, optional
            DESCRIPTION. The default is "exterior".
        c_0 : TYPE, optional
            DESCRIPTION. The default is 343.
        rho_0 : TYPE, optional
            DESCRIPTION. The default is 1.22.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        # inputs
        if isinstance(xSource, np.ndarray) is True:
            self.xSource = xSource
        else:
            self.xSource = np.array(xSource)
        self.xSystem = None
        self.QSource = np.array(volumeVelocity)
        self.QSystem = None
        self.frequency = frequency
        self.domain = domain
        self.c_0 = c_0
        self.rho_0 = rho_0
        self.kwargs = kwargs

        # check if mesh is v2
        self.meshPath = gtb.geometry.check_mesh(meshPath)
        if "radiatingElement" in kwargs:
            self.radiatingElement = kwargs["radiatingElement"]
        
        # other
        self.k = self.frequency * 2 * np.pi / self.c_0
        self.Ns = len(self.xSource)
        self.isComputed = False
        
        # initialize possible boundary conditions
        self.impedanceSurfaceIndex = []
        self.surfaceImpedance = []
        
        # get kwargs
        self.boundary_conditions = None
        self.direction = False
        self.tol = None
        self.parse_input()
        
        # initialize pressures and velocities arrays
        self.u_mesh = np.empty([len(frequency), self.Ns], dtype=object)  # separate sources
        self.u_mesh_Y = np.empty([len(frequency), self.Ns], dtype=object)
        self.p_mesh = np.empty([len(frequency), self.Ns], dtype=object)  # separate drivers
        self.u_total_mesh = np.empty([len(frequency)], dtype=object)   # summed sources
        self.p_total_mesh = np.empty([len(frequency)], dtype=object)   # summed sources
        self.propag_function = np.empty([len(frequency), self.Ns], dtype=object)
        
        # load simulation grid and mirror mesh if needed
        self.grid_sim = bempp.api.import_grid(self.meshPath)
        self.grid_init = bempp.api.import_grid(self.meshPath)
        self.grid_sim, self.sizeFactor = mirror_mesh(self.grid_init, self.boundary_conditions)
        self.vertices = np.shape(self.grid_sim.vertices)[1]
        
        if self.boundary_conditions is not None:
            self.mirror_system()
        else:
            self.xSystem = self.xSource
        
        
        # define space functions
        self.spaceP   = bempp.api.function_space(self.grid_sim, "DP", 0)
        self.identity = sparse.identity(self.spaceP, self.spaceP, self.spaceP)
        
        
        # Assign vibrometric coefficients to corresponding radiators // Assign surface velocity if not vibration data
        self.admittanceCoeff = getSurfaceAdmittance(self.impedanceSurfaceIndex, self.surfaceImpedance, 
                                                    frequency, self.spaceP, self.c_0, self.rho_0)        
        # driver reference
        self.LEM_enclosures = None
        self.radiator = None
        
        
    def parse_input(self):
        """
        Detects and assign kwargs to BEMOBJ variables.

        Returns
        -------
        None.

        """
        if "boundary_conditions" in self.kwargs:
            self.boundary_conditions = self.kwargs["boundary_conditions"].parameters
            self.initialize_conditions()
        elif "boundary_condition" in self.kwargs: # just in case user forget "s"
            self.boundary_conditions = self.kwargs["boundary_condition"].parameters
            self.initialize_conditions()
        if "direction" in self.kwargs:
            self.direction = self.kwargs["direction"]
        if "vibrometry_points" in self.kwargs:
            self.vibrometry_points = self.kwargs["vibrometry_points"]
        if "tol" in self.kwargs:
            self.tol = self.kwargs["tol"]
        else:
            self.tol = 1e-5
            
    def initialize_conditions(self):
        """
        Assign boundary conditions to the impedanceSurfaceIndex and surfaceImpedance 
        variables.

        Returns
        -------
        None.

        """
        
        for bc in self.boundary_conditions:
            if bc not in ["x", "X", "y", "Y", "z", "Z"]:
                self.impedanceSurfaceIndex.append(self.boundary_conditions[bc]["index"])
                self.surfaceImpedance.append(self.boundary_conditions[bc]["impedance"])
            else:
                pass
     
    
    def mirror_system(self):
        self.xSystem = np.copy(self.xSource)
        for param in self.boundary_conditions:
            if param in ["x", "X"]:
                xmirror = np.copy(self.xSystem)
                xmirror[:, 0] = 2*self.boundary_conditions.parameters[param]["offset"] - xmirror[:, 0]
                self.xSystem = np.concatenate([self.xSystem, xmirror])
            if param in ["y", "Y"]:
                ymirror = np.copy(self.xSystem)
                ymirror[:, 1] = 2*self.boundary_conditions.parameters[param]["offset"] - ymirror[:, 1]
                self.xSystem = np.concatenate([self.xSystem, ymirror])
            if param in ["z", "Z"]:
                zmirror = np.copy(self.xSystem)
                zmirror[:, 2] = 2*self.boundary_conditions.parameters[param]["offset"] - zmirror[:, 2]
                self.xSystem = np.concatenate([self.xSystem, zmirror])
        self.QSystem = np.tile(self.QSource, (self.sizeFactor, 1))
        
    
    def solve(self):
        """
        Compute the Boundary Element Method (BEM) solution for the loudspeaker system.

        This method performs the BEM computations to determine the acoustic pressure distribution on the mesh
        due to the contribution of individual speakers. The total pressure distribution is also computed by summing
        up the contributions of all speakers.

        Returns
        -------
        None

        Notes
        -----
        This function calculates the acoustic pressure distribution on the mesh at various frequencies using the
        Boundary Element Method. It iterates through the frequencies and speakers, calculates the necessary BEM
        operators (double layer, single layer), and uses GMRES solver to solve the BEM equation for the acoustic
        pressure distribution. The results are stored in class attributes for further analysis.

        """
        
        if self.domain == "exterior":
            domain_operator = -1
        elif self.domain == "interior":
            domain_operator = +1
        
        omega = 2 * np.pi * self.frequency
        k = -omega / self.c_0
        
        print("Computing pressure on mesh")
        if self.admittanceCoeff is None:
            for i in tqdm(range(len(k))):
                double_layer = helmholtz.double_layer(self.spaceP, self.spaceP,
                                                      self.spaceP, k[i])
                single_layer = helmholtz.single_layer(self.spaceP, self.spaceP,
                                                      self.spaceP, k[i])
                lhs = double_layer + domain_operator*(0.5 * self.identity)
                
                for rs in range(self.Ns):
                    xsource = self.xSystem[rs::self.Ns, :]
                    G  = greenMonopole(xsource, self.grid_sim, k[i])
                    n0 = incidenceCoeff(xsource, self.grid_sim, self.domain)
                    self.propag_function[i, rs] = buildGridFunction(self.spaceP, 
                                                        -1j*k[i]*n0*G)
                    rhs = single_layer * self.propag_function[i, rs]
                    self.u_mesh[i, rs], _ = gmres(lhs, rhs, tol=self.tol)
            self.isComputed = True
            
        elif self.admittanceCoeff is not None:
            for i in tqdm(range(len(k))):
                double_layer = helmholtz.double_layer(self.spaceP, self.spaceP,
                                                      self.spaceP, k[i])
                single_layer = helmholtz.single_layer(self.spaceP, self.spaceP,
                                                      self.spaceP, k[i])

                # ABSORBING SURFACES
                Yn = self.admittanceCoeff[:, i]  # all admittance coeff at current frequency
                yn_fun = bempp.api.GridFunction(self.spaceP, coefficients=Yn)  # ? doubts on its usefulness
                yn = DiagonalOperator(yn_fun.coefficients)
                
                lhs = ((double_layer + 
                       domain_operator*(0.5 * self.identity)).weak_form() - 
                       (1j*k[i]*single_layer.weak_form() * yn))
                
                for rs in range(self.Ns):
                    # PROPAGATION FROM SOURCE TO BOUNDARIES
                    xsource = self.xSystem[rs::self.Ns, :]
                    G  = greenMonopole(xsource, self.grid_sim, k[i])
                    n0 = incidenceCoeff(xsource, self.grid_sim, self.domain)
                    self.propag_function[i, rs] = buildGridFunction(self.spaceP, 
                                                        -1j*k[i]*n0*G)
                    
                    
                    rhs = (single_layer * self.propag_function[i, rs])
                    rhs = rhs.projections(self.spaceP)
                    
                    u_total_coeff, _ = scipy_gmres(lhs, rhs, rtol=self.tol)
                    u_total_coeff_Y = 1j*k[i]*u_total_coeff*yn_fun.coefficients
                    self.u_mesh[i, rs] = bempp.api.GridFunction(self.spaceP, 
                                                                coefficients=u_total_coeff)
                    self.u_mesh_Y[i, rs] = bempp.api.GridFunction(self.spaceP, 
                                                                coefficients=u_total_coeff_Y)

                    
            self.isComputed = True
        self.p_mesh = element2vertice_pressure(self.grid_sim, self.u_mesh)


    def getMicPressure(self, micPosition, individualSpeakers=False):
        """
        Get the pressure received at the considered microphones.

        Parameters
        ----------
        micPosition : numpy array
            Coordinates of the microphones (Cartesian). Shape: (nMic, 3)
        individualSpeakers : bool, optional
            If True, returns an array containing pressure received at each microphone from individual speakers.
            If False, returns the summed pressure received at each microphone from all speakers.
            Default is False.

        Returns
        -------
        pressure_mic : numpy array
            Pressure received at the specified microphones. Shape: (nFreq, nMic)

        Notes
        -----
        This function calculates the acoustic pressure received at the specified microphone positions for each frequency
        in the frequency array. The pressure is computed based on the BEM solution obtained from `computeBEM` method.
        It uses BEM operators (double layer, single layer) and their interactions with the speaker distributions to
        compute the pressure at each microphone position.

        """
        micPosition = np.array(micPosition).T
        nMic = np.shape(micPosition)[1]

        pressure_mic_array = np.zeros([len(self.frequency), nMic, self.Ns], dtype=complex)
        pressure_mic = np.zeros([len(self.frequency), nMic], dtype=complex)
        omega = 2 * np.pi * self.frequency
        k = -omega / self.c_0
        
        
        if self.admittanceCoeff is None:
            for i in tqdm(range(len(k))):
                for rs in range(self.Ns):
                    xsource = self.xSystem[rs::self.Ns, :]
                    single_pot = helmholtz_potential.single_layer(self.spaceP, 
                                                                  micPosition, k[i])
                    double_pot = helmholtz_potential.double_layer(self.spaceP, 
                                                                  micPosition, k[i])
                    for imag in range(len(xsource)):
                        dist = ((xsource[imag, 0]-micPosition[0, :])**2 + 
                                (xsource[imag, 1]-micPosition[1, :])**2 + 
                                (xsource[imag, 2]-micPosition[2, :])**2)**0.5
                        pressure_mic_array[i, :, rs] += np.reshape(np.exp(1j*k[i]*dist) / (4*np.pi*dist)
                             + double_pot.evaluate(self.u_mesh[i, rs])
                             - single_pot.evaluate(self.propag_function[i, rs]), nMic)
                    pressure_mic[i, :] += pressure_mic_array[i, :, rs]
            
            if individualSpeakers is True:
                out = (pressure_mic, pressure_mic_array)
            elif individualSpeakers is False:
                out = pressure_mic
        elif self.admittanceCoeff is not None:
            for i in tqdm(range(len(k))):
                for rs in range(self.Ns):
                    xsource = self.xSystem[rs::self.Ns, :]
                    single_pot = helmholtz_potential.single_layer(self.spaceP, 
                                                                  micPosition, k[i])
                    double_pot = helmholtz_potential.double_layer(self.spaceP, 
                                                                  micPosition, k[i])
                    for imag in range(len(xsource)):
                        dist = ((xsource[imag, 0]-micPosition[0, :])**2 + 
                                (xsource[imag, 1]-micPosition[1, :])**2 + 
                                (xsource[imag, 2]-micPosition[2, :])**2)**0.5
                        pressure_mic_array[i, :, rs] += np.reshape(np.exp(1j*k[i]*dist) / (4*np.pi*dist)
                             + double_pot.evaluate(self.u_mesh[i, rs])
                             - single_pot.evaluate(self.propag_function[i, rs])
                             + single_pot.evaluate(self.u_mesh_Y[i, rs]), nMic)
                    pressure_mic[i, :] += pressure_mic_array[i, :, rs]
            
            if individualSpeakers is True:
                out = (pressure_mic, pressure_mic_array)
            elif individualSpeakers is False:
                out = pressure_mic
        return out


#%% helper function
# def getSurfaceAdmittance(absorbingSurface, surfaceImpedance, freq, spaceP, c_0, rho_0):
#     """
#     Compute the total single layer coefficients linked to surfaces impedance.
#     :param absorbingSurface:
#     :param surfaceImpedance:
#     :param freq:
#     :param spaceP:
#     :return:
#     """
#     Nfft       = len(freq)
#     absSurf_in = np.array(absorbingSurface)
#     surfImp_in = surfaceImpedance
#     Nsurf      = len(absSurf_in)             # number of absorbing surfaces
#     grid       = spaceP.grid
#     dofCount   = spaceP.grid_dof_count

#     if absSurf_in.shape[0] == 0 :
#         admittanceMatrix = None
#     else:
#         admittanceMatrix = np.ones([dofCount, Nfft], dtype=complex) * 2.5e-3 # corresponds to 0.5% damping
#         for surf in range(Nsurf):
#             tmp_surf = absSurf_in[surf]  # current surface on which we apply admittance coefficients
#             vertex, _ = get_group_points(grid, tmp_surf)
#             for f in range(Nfft):
#                 try:
#                     Yn = rho_0 * c_0 / surfImp_in[surf][f] # rho_0 * c_0
#                 except:
#                     Yn = rho_0 * c_0 / surfImp_in[surf] # rho_0 * c_0
#                 admittanceMatrix[vertex, f] = np.ones(len(vertex)) * Yn
#     return admittanceMatrix      