#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:42:38 2023

@author: tom.munoz
"""
import os
import meshio
import bempp.api
from bempp.api.operators.boundary import helmholtz, sparse
from bempp.api.operators.potential import helmholtz as helmholtz_potential
from bempp.api.linalg import gmres
import numpy as np
import generalToolbox as gtb
import pyvista as pv
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore", message="splu requires CSC matrix format")
warnings.filterwarnings("ignore", message="splu converted its input to CSC format")
from electroacPy.globalVariables import air

# bempp.api.set_default_gpu_device_by_name('NVIDIA CUDA')
# bempp.api.BOUNDARY_OPERATOR_DEVICE_TYPE = 'gpu'
# bempp.api.POTENTIAL_OPERATOR_DEVICE_TYPE = 'gpu'
# bempp.api.DEFAULT_PRECISION = 'single'

class bem:
    def __init__(self, meshPath, radiatingSurface, surfaceVelocity,
                 freq_array, vibrometry_points=False, radiation_direction=False, boundary=False,
                 offset=False, backend='', c_0=air.c, rho_0=air.rho, tol=1e-5, domain="exterior"):
        """
        Create a bemppProblem object, using an input mesh, the index of the
        radiating surfaces (available in the mesh file) and corresponding
        velocities of these radiating surfaces.
        This class revolves around the simulation of loudspeaker systems.


        ! if "LogicError: create_buffer failed: INVALID_BUFFER_SIZE" appears,
        the surface given in radiatingSurface might not correspond to the one
        appearing in the mesh file !

        Parameters
        ----------
        meshPath : string
            path to the mesh file of the audio system to simulate.
        radiatingSurface : list,
            index of the radiating surfaces. ex:[2, 3, 4]
        surfaceVelocity : numpy array
            velocities of the radiating surfaces. shape: (nFreq, nSurface)
        freq_array: numpy array
            frequency axis corresponding to the velocities
        vibrometry_point: dict
            point cloud linked to vibrometric data. The default is False
        boundary : str, optional
            Indicate to which axis is the rigid boundary normal to. The default is False.
        offset : int or float, optional
            Boundary offset regarding its orientation axis. The default is False.
        backend : str, optional
            Choose in which software the radiation maps are shown. The default is ''.

        TODO
        - phase issues at low frequencies?
        """

        # identifier (save/load function)
        self.identifier = "BEM_EXT"
        self.rho_0 = rho_0
        self.c_0 = c_0
        self.tol = tol

        # check mesh version and display message if v4
        meshFile = open(meshPath)
        lines = meshFile.readlines()
        if lines[1][0] != '2':
            raise TypeError(
                "Mesh file is not in version 2. Errors will appear when mirroring mesh along boundaries.")
        meshFile.close()

        # check backend
        if backend != '':
            bempp.api.PLOT_BACKEND = backend

        # initialize data
        self.meshPath            = meshPath
        self.radiatingSurface    = np.array(radiatingSurface)
        self.surfaceVelocity     = surfaceVelocity
        self.nRad_S              = len(radiatingSurface)
        self.freq_array          = freq_array
        self.radiation_direction = radiation_direction
        self.isComputed          = False  # is the pressure computed on mesh

        # check surface velocity for vibrometric data
        self.isVibData = checkVelocityInput(surfaceVelocity)
        self.vibrometry_points = vibrometry_points

        # init pressures and velocities at zero
        self.u_mesh = np.empty([len(freq_array), self.nRad_S], dtype=object)  # separate drivers
        self.p_mesh = np.empty([len(freq_array), self.nRad_S], dtype=object)  # separate drivers
        self.u_total_mesh       = np.empty([len(freq_array)], dtype=object)   # summed drivers
        self.p_total_mesh       = np.empty([len(freq_array)], dtype=object)   # summed drivers

        # import grid and mirror it if floor_boundary != "none"
        self.grid      = bempp.api.import_grid(meshPath)
        self.grid_show = bempp.api.import_grid(meshPath)
        self.vertices  = np.shape(self.grid.vertices)[1]  # number of vertices
        self.boundary  = boundary
        self.offset    = offset

        # check for boundaries and mirror mesh if necessary
        grid_tot = bempp.api.import_grid(meshPath)
        gridPath = meshPath
        self.sizeFactor = 1

        # image / source
        if boundary is not False:
            for idx, val in enumerate(boundary):
                vertices = np.copy(grid_tot.vertices)
                elements = np.copy(grid_tot.elements)
                if val == "x" or val == "X":
                    if offset is not False:
                        vertices[0, :] = 2*offset[idx] - vertices[0, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                    else:
                        vertices[0, :] = vertices[0, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                elif val == "y" or val == "Y":
                    if offset is not False:
                        vertices[1, :] = 2*offset[idx] - vertices[1, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                    else:
                        vertices[1, :] = vertices[1, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                elif val == "z" or val == "Z":
                    if offset is not False:
                        vertices[2, :] = 2*offset[idx] - vertices[2, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                    else:
                        vertices[2, :] = vertices[2, :]
                        elements[[2, 0], :] = elements[[0, 2], :]
                else:
                    print("Instruction unclear, stuck the mesh up my RAM. Check 'x', 'y' or 'z' boundary definition(s).")
                grid_mirror = bempp.api.Grid(vertices, elements, domain_indices=grid_tot.domain_indices)
                grid_tot = bempp.api.grid.union([grid_tot, grid_mirror],
                                                [grid_tot.domain_indices, grid_mirror.domain_indices])
                self.sizeFactor *= 2
            self.grid = grid_tot
            self.vertices = np.shape(self.grid.vertices)[1]

        # define space functions
        self.spaceP   = bempp.api.function_space(self.grid, "P", 1)
        self.identity = sparse.identity(self.spaceP, self.spaceP, self.spaceP)

        # spaceU
        # DP, 0
        self.spaceU = bempp.api.function_space(self.grid, "DP", 0,
                                               segments=self.radiatingSurface)

        self.spaceU_list     = np.empty(self.nRad_S, dtype=object)
        self.u_callable_list = np.empty(self.nRad_S, dtype=object)
        
        # create a list of all velocity spaces (corresponding to radiators) as well as correction coefficients
        self.correctionCoefficients = []
        dof = np.zeros(self.nRad_S)
        for i in range(self.nRad_S):
            # create the spaceU function_space
            # DP, 0
            spaceU = bempp.api.function_space(self.grid, "DP", 0,
                                              segments=[radiatingSurface[i]])
            self.spaceU_list[i] = spaceU
            dof[i] = int(spaceU.grid_dof_count)   # degree of freedom of each radiators
            if isinstance(self.radiation_direction, bool) is False:  # if radiation_direction has been defined by user
                self.correctionCoefficients.append(getCorrectionCoefficients(spaceU.support_elements,
                                                             spaceU.grid.normals, radiation_direction[i]))
            else:
                self.correctionCoefficients.append(1)

        # Assign vibrometric coefficients to corresponding radiators // Assign surface velocity if not vibration data
        self.dof = dof
        maxDOF = int(np.max(dof) + 1)
        self.coeff_radSurf = np.zeros([len(freq_array), self.nRad_S, maxDOF], dtype=complex)
        for rs, isVib in enumerate(self.isVibData):
            if isVib is False:
                for f in range(len(self.freq_array)):
                    self.coeff_radSurf[f, rs, :int(dof[rs])] = surfaceVelocity[rs][f]
            if isVib is True:
                self.coeff_radSurf[:, rs, :int(dof[rs])] = getRadiationCoefficients(self.spaceU_list[rs].support_elements,
                                                                                    self.grid_show.centroids,
                                                                                    surfaceVelocity[rs],
                                                                                    vibrometry_points[rs],
                                                                                    self.sizeFactor)
        # driver reference
        self.LEM_enclosures = None
        self.radiator = None


    # %% MIRROR MESH
    def mirror_mesh(self, i, gridPath):
        """
        Mirror a mesh along a specified boundary plane.

        Parameters
        ----------
        i : int
            Index indicating the boundary plane to mirror along.
        gridPath : str
            Path to the input mesh file in Gmsh format.

        Returns
        -------
        mirrorPath : str
            Path to the output mirrored mesh file.

        Raises
        ------
        ValueError
            If the specified 'boundary' is not understood.

        Notes
        -----
        This function mirrors a given mesh along the specified boundary plane ('x', 'y', or 'z') by reversing the
        coordinates along the selected axis. It also reverses the normals of the cells to ensure correct orientation.

        """
        if self.boundary[i] == "x":
            ind_floor = 0
        elif self.boundary[i] == "y":
            ind_floor = 1
        elif self.boundary[i] == "z":
            ind_floor = 2
        else:
            raise ValueError("'boundary' not understood.")

        mesh = meshio.read(gridPath)
        points_n = mesh.points
        if self.offset is not False:
            points_n[:, ind_floor] = 2 * self.offset[i] - mesh.points[:, ind_floor]
        else:
            points_n[:, ind_floor] = - mesh.points[:, ind_floor]
        mesh.cells[0].data[:, [2, 0]] = mesh.cells[0].data[:, [0, 2]]  # will reverse normals
        mesh_mirror = meshio.Mesh(
            points_n,
            mesh.cells,
            # Optionally provide extra data on points, cells, etc.
            mesh.point_data,
            # Each item in cell data must match the cells array
            mesh.cell_data,
        )
        mirrorPath = "mirror_{}.msh".format(self.boundary)
        mesh_mirror.write(mirrorPath, binary=False, file_format="gmsh22")
        return mirrorPath


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
        omega = 2 * np.pi * self.freq_array
        k = -omega / self.c_0

        # individual speakers
        self.u_mesh = np.empty([len(k), self.nRad_S], dtype=object)
        self.p_mesh = np.empty([len(k), self.nRad_S], dtype=object)

        # sum of all speakers
        self.u_total_mesh = np.empty([len(k)], dtype=object)
        self.p_total_mesh = np.empty([len(k)], dtype=object)
        coefficients_summed = np.zeros([len(k), self.spaceP.grid_dof_count], dtype=complex)

        # error
        self.error = np.zeros([len(k), self.nRad_S])

        print("Computing pressure on mesh")
        for i in tqdm(range(len(k))):
            # creation of the double layer
            double_layer = helmholtz.double_layer(self.spaceP, self.spaceP,
                                                  self.spaceP, k[i])
            for rs in range(self.nRad_S):
                grid_dof = self.spaceU_list[rs].grid_dof_count
                
                # coeff_radSurf = np.ones(grid_dof, dtype=complex) * self.surfaceVelocity[i, rs]
                coeff_radSurf = self.coeff_radSurf[i, rs, :int(self.dof[rs])]

                # DP, 0
                spaceU = bempp.api.function_space(self.grid, "DP", 0,
                                                  segments=[self.radiatingSurface[rs]])

                # get velocity on current radiator
                u_total = bempp.api.GridFunction(spaceU, coefficients=-coeff_radSurf *
                                                                      self.correctionCoefficients[rs])

                # single layer
                single_layer = helmholtz.single_layer(spaceU,
                                                      self.spaceP, self.spaceP,
                                                      k[i])

                # pressure over the whole surface of the loudspeaker (p_total)
                p_total, _ = gmres(double_layer - 0.5 * self.identity,
                                   1j * omega[i] * self.rho_0 * single_layer * u_total,
                                   tol=self.tol, return_residuals=False)
                # self.error[i, rs] = e[-1]

                self.p_mesh[i, rs] = p_total # individual speakers
                self.u_mesh[i, rs] = u_total # individual speakers
                coefficients_summed[i, :] += p_total.coefficients
            self.p_total_mesh[i] = bempp.api.GridFunction(self.spaceP, coefficients=coefficients_summed[i])
        self.isComputed = True
        return None

    # %% COMPUTE PRESSURE AT MICROPHONE POSITIONS
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

        # check for CUDA
        # if "NVIDIA CUDA" in str(bempp.api.GPU_OPENCL_DRIVER_FOUND[1]):
        #     bempp.api.set_default_gpu_device_by_name("NVIDIA CUDA")
        #     bempp.api.POTENTIAL_OPERATOR_DEVICE_TYPE = "gpu"
        #     bempp.api.DEFAULT_PRECISION = 'single'
        #     print("NVIDIA CUDA")

        micPosition = np.array(micPosition).T
        nMic = np.shape(micPosition)[1]

        pressure_mic_array = np.zeros([len(self.freq_array), nMic, self.nRad_S], dtype=complex)
        pressure_mic = np.zeros([len(self.freq_array), nMic], dtype=complex)
        omega = 2 * np.pi * self.freq_array
        k = -omega / self.c_0


        print("\n" + "Computing pressure at microphones")
        for i in tqdm(range(len(k))):  # looping through frequencies
            for rs in range(self.nRad_S):
                # pressure received at microphones
                pressure_mic_array[i, :, rs] = np.reshape(
                    helmholtz_potential.double_layer(self.spaceP,
                                              micPosition,
                                              k[i]) * self.p_mesh[i, rs] \
                    - 1j * omega[i] * self.rho_0 * \
                    helmholtz_potential.single_layer(
                        self.spaceU_list[rs],
                        micPosition, k[i]) * self.u_mesh[i, rs], nMic)
                pressure_mic[i, :] += pressure_mic_array[i, :, rs]

        if individualSpeakers is True:
            out = (pressure_mic, pressure_mic_array)
        elif individualSpeakers is False:
            out = pressure_mic
        return out


    def get_surfacePressure(self, radiatingSurface):
        """
        Return the integration of the pressure over the corresponding surface
        """

        nodes_per_group = gtb.geometry.extract_bem_surface(self.grid)
        nodes_surface = nodes_per_group[str(radiatingSurface)]  # array
        pressure = np.zeros(len(self.freq_array), dtype=complex)
        for i in range(len(nodes_surface)):
            idx = nodes_surface[i]
            for f in range(len(self.freq_array)):
                idx_surf = np.argwhere(self.radiatingSurface == radiatingSurface)[0][0]
                pressure[f] += self.p_mesh[f][idx_surf].coefficients[idx]

        return pressure

    # def getMeshPressure(self, radiatingSurfaces="all", freq_idx=0):
    #     """
    #     Return pressure on mesh for selected surfaces
    #
    #     :param radiatingSurfaces:
    #     :return:
    #     """
    #     # print("PMESH:::: ", self.p_mesh[0])
    #     if radiatingSurfaces != "all":
    #         for i in range(len(radiatingSurfaces)):
    #             rad2add = np.argwhere(self.radiatingSurface == radiatingSurfaces[i])[0][0]
    #             print(rad2add)
    #
    #         p_out = np.zeros(self.spaceP.grid_dof_count, dtype=complex)
    #         for j in range(len(rad2add)):
    #             p_out += self.p_mesh[freq_idx][j].coefficients
    #     else:
    #         p_out = np.zeros(self.spaceP.grid_dof_count, dtype=complex)
    #         for i in range(self.nRad_S):
    #             p_out += self.p_mesh[freq_idx][i].coefficients
    #
    #
    #     return p_out


    ## %% DISPLAY VELOCITY ON MESH
    def displayVelocityMesh(self, frequency):
        """
        Display the velocity distribution on the mesh at a specific frequency.

        Parameters
        ----------
        frequency : float
            Frequency at which to display the velocity distribution.

        Returns
        -------
        plot_velocity : bempp.api.GridFunction
            GridFunction containing the velocity distribution on the mesh.

        Notes
        -----
        This function calculates and displays the velocity distribution on the mesh at a specific frequency.
        It uses the BEM++ library for visualization. The velocity values are displayed as real part.

        """
        indFreq, freqValue = gtb.findInArray(self.freq_array, frequency)
        plot_velocity = self.u_total_mesh[indFreq]
        return plot_velocity.plot(transformation='real')

    ##
    # %% REWORK OF THE DISPLAY PRESSURE AND DISPLAY MESH FUNCTIONS
    def displayMesh(self, edges=True, boundary=True):
        """
        Display the mesh of the system.

        Parameters
        ----------
        edges : bool, optional
            Whether to plot edges of the mesh. Default is True.
        boundary : bool, optional
            Whether to plot rigid boundaries. Default is True.

        Returns
        -------
        pv.Plotter
            A PyVista plotter displaying the mesh.

        Notes
        -----
        This function displays the mesh of the system using PyVista. It can visualize edges of the mesh and
        rigid boundaries if specified. Mesh vertices corresponding to different radiating surfaces are
        color-coded. Rigid boundaries (if defined) are displayed as light gray surfaces.

        """
        speakerMesh = pv.read(self.meshPath)

        # Map surface indices to colors
        speakerMesh.cell_colors = self.radiatingSurface

        pl = pv.Plotter(lighting='three lights')
        pl.add_mesh(speakerMesh, show_edges=edges, cmap='summer', show_scalar_bar=False, lighting=False)
        _ = pl.show_grid(color='k', fmt='%.3f')
        pl.add_axes(color='black')

        if self.boundary is not False:
            if boundary is True:
                for i in range(len(self.boundary)):
                    dimX = np.max(self.grid_show.vertices) * 2
                    dimY = np.max(self.grid_show.vertices) * 2
                    dimZ = np.max(self.grid_show.vertices) * 2
                    if self.boundary[i] == 'x':
                        vertices = np.array([[self.offset[i], -dimY, -dimZ],
                                             [self.offset[i], -dimY, dimZ],
                                             [self.offset[i], dimY, dimZ],
                                             [self.offset[i], dimY, -dimZ]])
                        face = [4, 0, 1, 2, 3]
                    elif self.boundary[i] == 'y':
                        vertices = np.array([[-dimX, self.offset[i], -dimZ],
                                             [-dimX, self.offset[i], dimZ],
                                             [dimX, self.offset[i], dimZ],
                                             [dimX, self.offset[i], -dimZ]])
                        face = [4, 0, 1, 2, 3]
                    elif self.boundary[i] == 'z':
                        vertices = np.array([[-dimX, -dimY, self.offset[i]],
                                             [-dimX, dimY, self.offset[i]],
                                             [dimX, dimY, self.offset[i]],
                                             [dimX, -dimY, self.offset[i]]])
                        face = [4, 0, 1, 2, 3]

                    surf = pv.PolyData(vertices, face)
                    pl.add_mesh(surf, color='lightgrey')
        # surf = pv.PolyData(self.grid.vertices.T, self.grid.elements)
        # pl.add_mesh(surf, color='lightgrey')
        pl.background_color = 'white'
        return pl.show()

    def displayPressureMesh(self):
        """
        Display the acoustic pressure distribution on the mesh using PyVista and interactive sliders.

        Returns
        -------
        pv.Plotter
            A PyVista plotter displaying the acoustic pressure distribution.

        Notes
        -----
        This function displays the acoustic pressure distribution on the mesh using PyVista. It creates an
        interactive slider widget allowing the user to visualize the pressure distribution at different
        frequencies. The SPL values are color-mapped and displayed on the mesh with an adjustable color bar.
        The frequency values for the slider are log-scaled.

        """
        pl = pv.Plotter(lighting='three lights')

        SPL_max_values = np.zeros(len(self.freq_array))
        for i in range(len(self.freq_array)):
            SPL_max_values = np.max(gtb.gain.SPL(self.p_total_mesh[i].coefficients))
        maxSPL = np.max(SPL_max_values) + 5
        minSPL = maxSPL - 45

        def create_mesh(value):
            # set log slider

            freq_value, frequencyValue = gtb.findInArray(self.freq_array, 10 ** (value))

            pl.remove_actor('SPL')

            sargs = dict(
                title_font_size=30,
                label_font_size=25,
                shadow=True,
                n_labels=5,
                italic=False,
                font_family="arial",
                color='k',
                title='SPL [dB] - ' + str(round(frequencyValue, 1)) + " Hz"
            )

            pl.background_color = 'lightyellow'

            speaker = pv.read(self.meshPath)
            pl.add_mesh(speaker, show_edges=True,
                        scalars=gtb.gain.SPL(self.p_total_mesh[int(freq_value)].coefficients), clim=[minSPL, maxSPL],
                        n_colors=20, show_scalar_bar=True, cmap="turbo",
                        scalar_bar_args=sargs, name='SPL', lighting=False)

            _ = pl.show_grid(color='k')
            pl.add_axes(color='black')

        slider = pl.add_slider_widget(
            create_mesh,
            rng=[np.log10(self.freq_array[0]), np.log10(self.freq_array[-1])],
            value=np.log10(250),
            title="Frequency [Hz]",
            # title_height=0.03,
            color='black',
            style='modern'
        )
        return pl.show()

    def plot_error(self):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        for rs in range(self.nRad_S):
            ax.semilogx(self.freq_array, self.error[:, rs], label='RS:{}'.format(rs))
        ax.grid(which="both")
        ax.set(xlabel="Frequency [Hz]", ylabel="Residuals")
        ax.legend(loc="best")
        plt.tight_layout()
        return None

## ==========================
# %% Pre-processing functions

def checkVelocityInput(surfaceVelocity):
    """
    Check the surfaceVelocity input parameter for vibrometric data: if surfaceVelocity[i] is 1 dimensional: velocity,
    if surfaceVelocity[i] is 2 dimensional: vibrometric data.
    :param surfaceVelocity:
    :return:
    """
    # velocity_out   = []
    # vibrometry_out = []
    isVibData      = []
    for i in range(len(surfaceVelocity)):
        if len(surfaceVelocity[i].shape) == 1:
            # velocity_out.append(surfaceVelocity[i])
            isVibData.append(False)
        elif len(surfaceVelocity[i].shape) == 2:
            # vibrometry_out.append(surfaceVelocity[i])
            isVibData.append(True)

    # velocity_out   = np.array(velocity_out).T
    return isVibData

def getRadiationCoefficients(support_elements, centroids, vibrometric_data, vibrometry_points, sizeFactor):
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

    if sizeFactor > 1:
        coefficients = np.tile(coefficients, (1, sizeFactor))
        # print(coefficients.shape)
    return coefficients

def getCorrectionCoefficients(support_elements, normals, radiation_direction):
    """
    Depending on the radiator's defined radiation direction, this function returns coefficients that are to be applied
    on the triangles of the considered radiator. The user can set specific direction on a single radiator only: by
    setting corresponding direction to False (direction will thus be normal to elementary surfaces)
    :param support_elements:
    :param normals:
    :param radiation_direction:
    :return:
    """
    correctionCoefficients = np.zeros(len(support_elements))

    if isinstance(radiation_direction, bool) is False: # check for cases like [[1, 0, 0], False, [1, 0, 0]]
        for i, element in enumerate(support_elements):
            correctionCoefficients[i] = (np.dot(radiation_direction, normals[element, :]) /
                                               (np.linalg.norm(normals[element, :]) * np.linalg.norm(radiation_direction)))
    else: # if [x, y, z] norm is not defined, set radiation to be normal to elementary surfaces
        correctionCoefficients = 1
    return correctionCoefficients


## coefficients / elements attribution
#     def assignUFFdata(self, uffFile, rotation=[0, 0, 0]):
#         freq, acc, R = gtb.io.loadUFF(uffFile)
#         self.freq_array, idx_d = gtb.freqop.decimate_frequency_axis(freq, self.freq_array)
#         acc_d = acc[:, idx_d]
#         R = gtb.geometry.rotatePointCloud(R, rotation_angles_deg=rotation)
#         uff_coefficients = []
#         for i in self.radiatingSurface:
#             index_rad_surface = self.spaceU_list[i].localised_space.support_elements
#             vertice_location = self.grid.centroids[list(index_rad_surface), :]
#             vertice_location_center = gtb.geometry.recenterZero(vertice_location)
#             uff_coefficients.append(getCoeff2Points(centroids_mesh,
#                                               centroids_measure,
#                                               acc_d, self.freq_array))


# ## Helpful functions
# def getCoeff2Points(centroids_mesh, centroids_measure, acc_data, freq_data):
#     ind = np.zeros(len(centroids_mesh))
#     val = np.zeros([len(centroids_mesh), 3])
#     ind_vert = []
#     for i in range(len(centroids_mesh)):
#         r = centroids_mesh[i, :]
#         index, value = gtb.geometry.findClosestPoint(centroids_measure, r[0], r[1], r[2])
#         ind[i] = index  # Rc[index] that corresponds to the closest point to vertice_location_center
#         val[i, :] = value
#     coefficients = np.zeros([len(centroids_mesh), Nfft], dtype=complex)
#     for i in range(len(centroids_mesh)):
#         numPoints = len(np.argwhere(ind == ind[i]))
#         coefficients[i, :] = acc_data[int(ind[i]), :] / gtb.freqop.laplace(freq_data) / numPoints
#     return coefficients