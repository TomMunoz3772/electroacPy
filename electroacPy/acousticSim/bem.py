#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:42:38 2023

@author: tom.munoz
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
# from pyopencl import CompilerWarning
import electroacPy.general as gtb
from .ACSHelpers_ import getSurfaceAdmittance

warnings.filterwarnings("ignore", message="splu requires CSC matrix format")
warnings.filterwarnings("ignore", message="splu converted its input to CSC format")
# warnings.filterwarnings("ignore", category=CompilerWarning)

# bempp.api.set_default_gpu_device_by_name('NVIDIA CUDA')
# bempp.api.BOUNDARY_OPERATOR_DEVICE_TYPE = 'gpu'
# bempp.api.POTENTIAL_OPERATOR_DEVICE_TYPE = 'gpu'
# bempp.api.DEFAULT_PRECISION = 'single'

try:
    from pyopencl import CompilerWarning
    warnings.filterwarnings("ignore", category=CompilerWarning)
except:
    None


class bem:
    def __init__(self, meshPath, radiatingElement, velocity, frequency, 
                 domain="exterior", c_0=343, rho_0=1.22, **kwargs):
        """
        Create a BEM object.

        Parameters
        ----------
        meshPath : str
            Path to simulation mesh.
        radiatingElement : list of int
            Reference to BEM physical group(s) - should be surfaces.
        velocity : list of numpy array
            Velocity of radiating surfaces.
        frequency : numpy array
            Range of simulation.
        domain : str, optional
            Domain of simulation: "interior" or "exterior". 
            Will change the BEM equation. Default is "exterior".
        c_0 : float, optional
            Speed of sound in the propagation medium. Default is c=343 m/s (air).
        rho_0 : float, optional
            Medium density. Default is rho=1.22 kg/m^3 (air).
        
        **kwargs
        --------
        boundary_condition(s) : boundaryCondition object,
            Define boundary conditions on BEM mesh (infinite 
            reflection plane, absorption surface, etc.)
        direction : list directional vectors, 
            Use this kwarg to set non-normal velocity on mesh. ex: [[1, 0, 0], [1, 0, 0]]
            for two radiating surfaces with velocity toward +x.
        vibrometry_points : ndarray,
            Position of measured vibrometry points. Should be of size (Npoints, 3).                
        tol : float
            Tolerance of the GMRES solver. By default is 1E-5.

        Returns
        -------
        None.

        """
        # get main parameters
        self.radiatingElement = radiatingElement
        self.velocity = velocity
        self.frequency = frequency
        self.c_0 = c_0
        self.rho_0 = rho_0
        self.domain = domain
        self.kwargs = kwargs
        
        # check if mesh is v2
        self.meshPath = gtb.geometry.check_mesh(meshPath)
        
        # initialize possible boundary conditions
        self.impedanceSurfaceIndex = []
        self.surfaceImpedance = []
        
        # get kwargs
        self.boundary_conditions = None
        self.direction = False
        self.vibrometry_points = None
        self.tol = None
        self.parse_input()
        
        # other parameters
        self.Ns = len(self.radiatingElement)
        self.isComputed = False
        self.is2dData = checkVelocityInput(self.velocity)
        
        # initialize pressures and velocities arrays
        self.u_mesh = np.empty([len(frequency), self.Ns], dtype=object)  # separate sources
        self.p_mesh = np.empty([len(frequency), self.Ns], dtype=object)  # separate drivers
        self.u_total_mesh = np.empty([len(frequency)], dtype=object)   # summed sources
        self.p_total_mesh = np.empty([len(frequency)], dtype=object)   # summed sources
        
        # load simulation grid and mirror mesh if needed
        self.grid_sim = bempp.api.import_grid(self.meshPath)
        self.grid_init = bempp.api.import_grid(self.meshPath)
        self.grid_sim, self.sizeFactor = mirror_mesh(self.grid_init, 
                                                     self.boundary_conditions)
        self.vertices = np.shape(self.grid_sim.vertices)[1]

        # define space functions
        self.spaceP   = bempp.api.function_space(self.grid_sim, "P", 1)
        self.identity = sparse.identity(self.spaceP, self.spaceP, self.spaceP)

        # initialize flow space
        self.spaceU_freq     = np.empty(self.Ns, dtype=object)
        self.u_callable_freq = np.empty(self.Ns, dtype=object)
        
        # create a list of all velocity spaces (corresponding to radiators) as well as correction coefficients
        self.correctionCoefficients = []
        dof = np.zeros(self.Ns)
        for i in range(self.Ns):
            spaceU = bempp.api.function_space(self.grid_sim, "DP", 0,
                                              segments=[radiatingElement[i]])
            self.spaceU_freq[i] = spaceU
            dof[i] = int(spaceU.grid_dof_count)   # degree of freedom of each radiators
            if isinstance(self.direction, bool) is False:  # if direction has been defined by user
                self.correctionCoefficients.append(getCorrectionCoefficients(spaceU.support_elements,
                                                   spaceU.grid.normals, self.direction[i]))
            else:
                self.correctionCoefficients.append(1)
        
        # Assign vibrometric coefficients to corresponding radiators // Assign surface velocity if not vibration data
        self.dof = dof
        maxDOF = int(np.max(dof) + 1)
        self.coeff_radSurf = np.zeros([len(frequency), self.Ns, maxDOF], dtype=complex)
        for rs, isVib in enumerate(self.is2dData):
            if isVib is False: # assign velocity
                for f in range(len(self.frequency)): # TODO: maybe try to remove that loop
                    self.coeff_radSurf[f, rs, :int(dof[rs])] = velocity[rs][f]
            if isVib is True:  # assign velocity from vibration 
                self.coeff_radSurf[:, rs, :int(dof[rs])] = getRadiationCoefficients(self.spaceU_freq[rs].support_elements,
                                                                                    self.grid_init.centroids,
                                                                                    velocity[rs],
                                                                                    self.vibrometry_points[rs],
                                                                                    self.sizeFactor)
        
        self.admittanceCoeff = getSurfaceAdmittance(self.impedanceSurfaceIndex, self.surfaceImpedance, 
                                                    frequency, self.spaceP, self.c_0, self.rho_0)        
        # driver reference
        self.LEM_enclosures = None
        self.radiator = None

    def discard_frequency(self):
        return None
        
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

        # individual speakers
        self.u_mesh = np.empty([len(k), self.Ns], dtype=object)
        self.p_mesh = np.empty([len(k), self.Ns], dtype=object)

        # sum of all speakers
        self.u_total_mesh = np.empty([len(k)], dtype=object)
        self.p_total_mesh = np.empty([len(k)], dtype=object)

        # error
        self.error = np.zeros([len(k), self.Ns])

        print("Computing pressure on mesh")
        if self.admittanceCoeff is None:
            for i in tqdm(range(len(k))):
                # creation of the double layer
                double_layer = helmholtz.double_layer(self.spaceP, self.spaceP,
                                                      self.spaceP, k[i])
                lhs = double_layer + 0.5 * self.identity * domain_operator
                for rs in range(self.Ns):                
                    coeff_radSurf = self.coeff_radSurf[i, rs, :int(self.dof[rs])]
                    spaceU = bempp.api.function_space(self.grid_sim, "DP", 0,
                                                      segments=[self.radiatingElement[rs]])
    
                    # get velocity on current radiator
                    u_total = bempp.api.GridFunction(spaceU, 
                                                     coefficients=-coeff_radSurf *
                                                     self.correctionCoefficients[rs])
                    # single layer
                    single_layer = helmholtz.single_layer(spaceU,
                                                          self.spaceP, self.spaceP,
                                                          k[i])
    
                    # pressure over the whole surface of the loudspeaker (p_total)
                    rhs = 1j * omega[i] * self.rho_0 * single_layer * u_total
                    p_total, _ = gmres(lhs, rhs, tol=self.tol, 
                                       return_residuals=False)
    
                    self.p_mesh[i, rs] = p_total # individual speakers
                    self.u_mesh[i, rs] = u_total # individual speakers
            self.isComputed = True
            
        elif self.admittanceCoeff is not None:
            for i in tqdm(range(len(k))):
                # creation of the double layer
                double_layer = helmholtz.double_layer(self.spaceP, self.spaceP,
                                                      self.spaceP, k[i])
                # admittance single layer
                single_layer_Y = helmholtz.single_layer(self.spaceP, self.spaceP,
                                                        self.spaceP, k[i])
                # ABSORBING SURFACES
                Yn = self.admittanceCoeff[:, i]  # all admittance coeff at current frequency
                yn_fun = bempp.api.GridFunction(self.spaceP, coefficients=Yn)  # ? doubts on its usefulness
                yn = DiagonalOperator(yn_fun.coefficients)

                # building left-hand-side term
                lhs = ((double_layer + 0.5*self.identity * domain_operator).weak_form()
                       - (1j*k[i]*single_layer_Y.weak_form()*yn))    
                
                for rs in range(self.Ns):
                    # RADIATING SURFACES
                    coeff_radSurf = self.coeff_radSurf[i, rs, :int(self.dof[rs])]

                    spaceU = bempp.api.function_space(self.grid_sim, "DP", 0,
                                                      segments=[self.radiatingElement[rs]])

                    # get velocity on current radiator
                    u_total = bempp.api.GridFunction(spaceU, 
                                                     coefficients=-coeff_radSurf *
                                                     self.correctionCoefficients[rs])

                    # single layer - radiating surface
                    single_layer = helmholtz.single_layer(spaceU,
                                                          self.spaceP, self.spaceP,
                                                          k[i])


                    rhs = 1j * omega[i] * self.rho_0 * single_layer * u_total
                    rhs = rhs.projections(self.spaceP)

                    # pressure over the whole surface of the loudspeaker (p_total)
                    p_total_coefficients, _ = scipy_gmres(lhs, rhs, rtol=self.tol)
                    p_total = bempp.api.GridFunction(self.spaceP, coefficients=p_total_coefficients)

                    self.p_mesh[i, rs] = p_total  # individual speakers
                    self.u_mesh[i, rs] = u_total  # individual speakers
            self.isComputed = True 
        return None
        
    
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

        print("\n" + "Computing pressure at microphones")
        for i in tqdm(range(len(k))):  # looping through frequencies
            for rs in range(self.Ns):
                # pressure received at microphones
                pressure_mic_array[i, :, rs] = np.reshape(
                    helmholtz_potential.double_layer(self.spaceP,
                                              micPosition,
                                              k[i]) * self.p_mesh[i, rs] \
                    - 1j * omega[i] * self.rho_0 * \
                    helmholtz_potential.single_layer(
                        self.spaceU_freq[rs],
                        micPosition, k[i]) * self.u_mesh[i, rs], nMic)
                pressure_mic[i, :] += pressure_mic_array[i, :, rs]

        if individualSpeakers is True:
            out = (pressure_mic, pressure_mic_array)
        elif individualSpeakers is False:
            out = pressure_mic
        return out



# %%useful functions
def checkVelocityInput(velocity):
    """
    Check the velocity input parameter for vibrometric data: if velocity[i] is 1 dimensional: velocity,
    if velocity[i] is 2 dimensional: vibrometric data.
    :param velocity:
    :return:
    """
    isVibData      = []
    for i in range(len(velocity)):
        if len(velocity[i].shape) == 1:
            isVibData.append(False)
        elif len(velocity[i].shape) == 2:
            isVibData.append(True)
    return isVibData


def mirror_mesh(grid_init, boundary_conditions):
    bc = boundary_conditions
    size_factor = 1
    grid_tot = grid_init
    if bc is not None:
        for item in bc:
            boundary = bc[item]
            if boundary["type"] == "infinite_baffle":
                offset = boundary["offset"]
                vertices = np.copy(grid_tot.vertices)
                elements = np.copy(grid_tot.elements)
                if item in ["x", "X"]:
                    if offset != 0:
                        vertices[0, :] = 2*offset - vertices[0, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                    else:
                        vertices[0, :] = vertices[0, :]
                        elements[[2, 0], :] = elements[[0, 2], :]                    
                elif item in ["y", "Y"]:
                    if offset != 0 :
                        vertices[1, :] = 2*offset - vertices[1, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                    else:
                        vertices[1, :] = vertices[1, :]
                        elements[[2, 0], :] = elements[[0, 2], :]    
                elif item in ["z", "Z"]:
                    if offset is not False:
                        vertices[2, :] = 2*offset - vertices[2, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                    else:
                        vertices[2, :] = vertices[2, :]
                        elements[[2, 0], :] = elements[[0, 2], :]               
                else:
                    print("{} not infinite baffle.".format(item))
                grid_mirror = bempp.api.Grid(vertices, elements, domain_indices=grid_tot.domain_indices)
                grid_tot = bempp.api.grid.union([grid_tot, grid_mirror],
                                                [grid_tot.domain_indices, grid_mirror.domain_indices])
                size_factor *= 2
    return grid_tot, size_factor


#%% Radiation coefficients
def getRadiationCoefficients(support_elements, centroids, vibrometric_data, 
                             vibrometry_points, sizeFactor):
    """
    Build radiation coefficients from vibrometric dataset.

    :param support_elements:
    :param centroids:
    :param vibrometric_data:
    :return:
    """
    
    if sizeFactor > 1:
        slices = gtb.slice_array_into_parts(support_elements, sizeFactor)
        vertex_location = centroids[slices[0], :]          # get [x, y, z] vertex position of radiator
        vertex_center   = gtb.geometry.recenterZero(vertex_location)    # recenter geometry at [x=0, y=0, z=0]
    else:
        vertex_location = centroids[support_elements, :]  # get [x, y, z] vertex position of radiator
        vertex_center = gtb.geometry.recenterZero(vertex_location)  # recenter geometry at [x=0, y=0, z=0]

    Nfft = np.shape(vibrometric_data)[1]
    _, coefficients = gtb.geometry.points_within_radius(vertex_center,
                                                        vibrometry_points,
                                                        5e-3, vibrometric_data, Nfft)

    print("COEFFICIENTS: ", coefficients.shape)
    if sizeFactor > 1:
        coefficients = np.tile(coefficients, (1, sizeFactor))
    return coefficients

def getCorrectionCoefficients(support_elements, normals, direction):
    """
    Depending on the radiator's defined radiation direction, this function returns coefficients that are to be applied
    on the triangles of the considered radiator. The user can set specific direction on a single radiator only: by
    setting corresponding direction to False (direction will thus be normal to elementary surfaces)
    :param support_elements:
    :param normals:
    :param direction:
    :return:
    """
    correctionCoefficients = np.zeros(len(support_elements))

    if isinstance(direction, bool) is False: # check for cases like [[1, 0, 0], False, [1, 0, 0]]
        for i, element in enumerate(support_elements):
            correctionCoefficients[i] = (np.dot(direction, normals[element, :]) /
                                               (np.linalg.norm(normals[element, :]) * np.linalg.norm(direction)))
    else: # if [x, y, z] norm is not defined, set radiation to be normal to elementary surfaces
        correctionCoefficients = 1
    return correctionCoefficients


#%% Impedance functions
def get_group_points(grid, group_number):
    domain_indices = grid.domain_indices
    elements       = grid.elements
    
    # Find the indices where domain_indices match the group_number
    group_indices = np.where(domain_indices == group_number)  # indices of support elements

    # Get the corresponding columns of elements
    group_points = elements[:, group_indices]  # segments (and NOT points!) part of group_number

    # Reshape to a 1D array and use unique to get unique values
    unique_points   = np.unique(group_points)

    return unique_points, group_indices

#%% boundary conditions
class boundaryConditions:
    def __init__(self, rho=1.22, c=343):
        self.parameters = {}
        self.c = c
        self.rho = rho
        self.Zc = rho*c
    
    def addInfiniteBoundary(self, normal, offset=0, **kwargs):
        if normal not in ["x", "y", "z", "X", "Y", "Z"]:
            raise ValueError("Normal to axis should be x, y or z.")
        self.parameters[normal] = {}
        self.parameters[normal]["offset"] = offset
        self.parameters[normal]["type"] = "infinite_baffle"
        if "absorption" in kwargs:
            self.parameters[normal]["absorption"] = kwargs["absorption"]
        elif "impedance" in kwargs:
            self.parameters[normal]["impedance"] = kwargs["impedance"]
        if "frequency" in kwargs:
            self.parameters[normal]["frequency"] = kwargs["frequency"]
        else:
            None
            
    def addSurfaceImpedance(self, name, index, data_type, value,
                            frequency=None, targetFrequency=None,
                            interpolation="linear"):
        """
        Add an impedance to a surface.

        Parameters
        ----------
        name : str,
            name of surface.
        index : int
            reference to BEM mesh.
        data_type : float or array of float
            Type of data: "impedance", "absorption", "reflection", "admittance".
        value : numpy array
            Data to apply on surface.
        frequency : numpy array, optional
            Range of impedance data (if using frequency-dependant impedance).
            The default is None.
        targetFrequency : numpy array, optional
            BEM frequency range - use it if impedance data as a different frequency-axis
            than the BEM study. The default is None.
        interpolation : str, optional
            How to interpolate impedance data on frequency range. Either "linear" or "cubic".
            The default is "linear".
            
        Returns
        -------
        None.
        
        Notes
        -----
        If data_type is "impedance", the corresponding value should be given without 
        normalization to characteristic impedance of air.

        """
        
        self.parameters[name] = {}
        self.parameters[name]["index"] = index
        self.parameters[name]["data_type"] = data_type
        self.parameters[name]["type"] = "absorption_surface"
        self.parameters[name]["frequency"] = frequency
        
        if data_type == "impedance":
            Z = self.parameters[name]["impedance"]
            self.parameters[name]["impedance"] = value
            self.parameters[name]["admittance"] = 1/Z
            self.parameters[name]["absorption"] = 1 - np.abs((Z-self.Zc)/(Z+self.Zc))**2
        elif data_type == "reflection":
            self.parameters[name]["reflection"] = value
            self.parameters[name]["impedance"] = (1+value) / (1-value) * self.Zc
            self.parameters[name]["admittance"] = 1/self.parameters[name]["impedance"]
        elif data_type == "absorption":
            self.parameters[name]["absorption"] = value
            self.parameters[name]["impedance"] = (2-value) / value *  \
                                                self.Zc
            self.parameters[name]["admittance"] = 1/self.parameters[name]["impedance"]
        elif data_type == "admittance":
            self.parameters[name]["admittance"] = value
            self.parameters[name]["impedance"] = 1 / value
        else:
            raise ValueError("'data_type' not understood. Try 'impedance', " +
                             "'admittance', 'reflection' or 'absorption'")
        
        if isinstance(frequency, np.ndarray) is True and isinstance(value, np.ndarray) is True:
            if len(frequency) != len(value):
                raise Exception("'value' and 'frequency' should have the same size.")
            
            if targetFrequency is None:
                self.parameters[name]["frequency"] = frequency
            else:
                if interpolation in ["linear", None]:
                    self.parameters[name]["impedance"] = np.interp(targetFrequency, 
                                                                   frequency, 
                                                                   self.parameters[name]["impedance"])
                    self.parameters[name]["admittance"] = np.interp(targetFrequency, 
                                                                    frequency, 
                                                                    self.parameters[name]["admittance"])
                    self.parameters[name]["frequency"] = targetFrequency
                elif interpolation == "cubic":
                    from scipy.interpolate import CubicSpline
                    csImp = CubicSpline(frequency, self.parameters[name]["impedance"])
                    csAdm = CubicSpline(frequency, self.parameters[name]["admittance"])
                    self.parameters[name]["impedance"] = csImp(targetFrequency)
                    self.parameters[name]["admittance"] = csAdm(targetFrequency)
                    self.parameters[name]["frequency"] = targetFrequency
                