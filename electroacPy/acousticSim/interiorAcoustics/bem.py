import os
import meshio
import bempp.api
from bempp.api.operators.boundary import helmholtz, sparse
from bempp.api.operators.potential import helmholtz as helmholtz_potential
from bempp.api.assembly.discrete_boundary_operator import DiagonalOperator
from bempp.api.linalg import gmres
from scipy.sparse.linalg import gmres as scipy_gmres
import numpy as np
import generalToolbox as gtb
import pyvista as pv
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore", message="splu requires CSC matrix format")
warnings.filterwarnings("ignore", message="splu converted its input to CSC format")
from electroacPy.globalVariables import air


class bem:
    def __init__(self, meshPath, radiatingSurface, surfaceVelocity,
                 freq_array, absorbingSurface=[], surfaceImpedance=[],
                 vibrometry_points=False, radiation_direction=False, backend='',
                 c_0=air.c, rho_0=air.rho):
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
        boundary : str, optional
            Indicate to which axis is the rigid boundary normal to. The default is "none".
        offset : int or float, optional
            Boundary offset regarding its orientation axis. The default is "none".
        backend : str, optional
            Choose in which software the radiation maps are shown. The default is ''.

        """

        # identifier (save/load function)
        self.identifier = "BEM_INT"
        self.c_0 = c_0
        self.rho_0 = rho_0

        # check mesh version and display message if v4
        meshFile = open(meshPath)
        lines = meshFile.readlines()
        if lines[1][0] != '2':
            raise TypeError(
                "Mesh file is not in version 2. Errors will appear when mirroring mesh along boundaries.")
        meshFile.close()

        if backend != '':
            bempp.api.PLOT_BACKEND = backend
        self.meshPath = meshPath
        self.radiatingSurface = np.array(radiatingSurface)
        self.surfaceVelocity = surfaceVelocity #np.array(surfaceVelocity).T
        self.nRad_S = len(radiatingSurface)
        self.absorbingSurface = absorbingSurface
        self.surfaceImpedance = surfaceImpedance
        self.radiation_direction = radiation_direction
        self.freq_array = freq_array
        self.isComputed = False  # is the pressure computed on mesh

        # check surface velocity for vibrometric data
        self.isVibData = checkVelocityInput(surfaceVelocity)
        self.vibrometry_points = vibrometry_points

        # init pressures and velocities at zero
        self.u_mesh = np.empty([len(freq_array), self.nRad_S], dtype=object)
        self.p_mesh = np.empty([len(freq_array), self.nRad_S], dtype=object)
        self.u_total_mesh = np.empty([len(freq_array)], dtype=object)
        self.p_total_mesh = np.empty([len(freq_array)], dtype=object)

        # init acoustical impedance to zero
        # all driver applied for each
        self.Zab_total = np.zeros([len(freq_array), len(radiatingSurface)], dtype=complex)
        # separated impedance (single speaker radiating in enclosure)
        self.Zab_single = np.zeros([len(freq_array), len(radiatingSurface)], dtype=complex)

        # import grid and mirror it if floor_boundary != "none"
        self.grid = bempp.api.import_grid(meshPath)
        self.grid_show = bempp.api.import_grid(meshPath)
        self.vertices = np.shape(self.grid.vertices)[1]  # number of vertices

        # define some space functions
        self.spaceP = bempp.api.function_space(self.grid, "P", 1)
        self.identity = sparse.identity(self.spaceP, self.spaceP, self.spaceP)


        # INPUT VELOCITY AND SOURCE'S RADIATION COEFFICIENTS
        # create a list of all velocity spaces (corresponding to radiators) as well as correction coefficients
        self.spaceU_list = np.empty(self.nRad_S, dtype=object)
        self.correctionCoefficients = []
        dof = np.zeros(self.nRad_S)
        for i in range(self.nRad_S):
            # create the spaceU function_space
            spaceU = bempp.api.function_space(self.grid, "DP", 0,
                                              segments=[radiatingSurface[i]])
            self.spaceU_list[i] = spaceU
            dof[i] = int(spaceU.grid_dof_count)  # degree of freedom of each radiators
            if isinstance(self.radiation_direction, bool) is False:  # if radiation_direction has been defined by user
                self.correctionCoefficients.append(getCorrectionCoefficients(spaceU.support_elements,
                                                                             spaceU.grid.normals,
                                                                             radiation_direction[i]))
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
                self.coeff_radSurf[:, rs, :int(dof[rs])] = getRadiationCoefficients(
                    self.spaceU_list[rs].support_elements,
                    self.grid_show.centroids,
                    surfaceVelocity[rs],
                    vibrometry_points[rs], 1)


        # ADMITTANCE OF ABSORBING SURFACES - SINGLE LAYER COMPUTATIONS
        # spaceY - one total singleLayer per frequency
        self.admittanceCoeff = getSurfaceAdmittance(absorbingSurface, surfaceImpedance, freq_array,
                                                    self.spaceP, self.c_0, self.rho_0)

        # driver reference
        self.LEM_enclosures = None
        self.radiator = None

        # link to other studies
        self.exterior_link = False

        # "False" variables for save function



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
        k = - omega / self.c_0

        # individual speakers
        self.u_mesh = np.empty([len(k), self.nRad_S], dtype=object)
        self.p_mesh = np.empty([len(k), self.nRad_S], dtype=object)

        # sum of all speakers
        self.u_total_mesh = np.empty([len(k)], dtype=object)
        self.p_total_mesh = np.empty([len(k)], dtype=object)
        coefficients_summed = np.zeros([len(k), self.spaceP.grid_dof_count], dtype=complex)

        print("Computing pressure on mesh")
        for i in tqdm(range(len(k))):
            # creation of the double layer
            double_layer = helmholtz.double_layer(self.spaceP, self.spaceP,
                                                  self.spaceP, k[i])
            # admittance single layer
            single_layer_Y = helmholtz.single_layer(self.spaceP, self.spaceP,
                                                    self.spaceP, k[i])
            for rs in range(self.nRad_S):
                # grid_dof = self.spaceU_list[rs].grid_dof_count

                # RADIATING SURFACES
                coeff_radSurf = self.coeff_radSurf[i, rs, :int(self.dof[rs])]

                spaceU = bempp.api.function_space(self.grid, "DP", 0,
                                                  segments=[self.radiatingSurface[rs]])

                # get velocity on current radiator
                u_total = bempp.api.GridFunction(spaceU, coefficients=-coeff_radSurf *
                                                                      self.correctionCoefficients[rs])

                # single layer - radiating surface
                single_layer = helmholtz.single_layer(spaceU,
                                                      self.spaceP, self.spaceP,
                                                      k[i])

                # ABSORBING SURFACES
                Yn = self.admittanceCoeff[:, i]  # all admittance coeff at current frequency
                yn_fun = bempp.api.GridFunction(self.spaceP, coefficients=Yn)  # ? doubts on its usefulness
                yn = DiagonalOperator(yn_fun.coefficients)

                # building equations
                lhs = (double_layer + 0.5*self.identity).weak_form() - (1j*k[i]*single_layer_Y.weak_form()*yn)
                rhs = 1j * omega[i] * self.rho_0 * single_layer * u_total
                rhs = rhs.projections(self.spaceP)

                # pressure over the whole surface of the loudspeaker (p_total)
                p_total_coefficients, _ = scipy_gmres(lhs, rhs, rtol=1E-5)
                p_total = bempp.api.GridFunction(self.spaceP, coefficients=p_total_coefficients)

                self.p_mesh[i, rs] = p_total  # individual speakers
                self.u_mesh[i, rs] = u_total  # individual speakers
                coefficients_summed[i, :] += p_total.coefficients
            self.p_total_mesh[i] = bempp.api.GridFunction(self.spaceP, coefficients=coefficients_summed[i])

        # self.getAcousticalImpedance()
        self.isComputed = True
        return None


    # Solve acoustical impedance on radiating surfaces
    def getAcousticalImpedance(self):
        """
        Compute the acoustical impedance applied on the radiating surfaces of the BEM object.

        :return:
        None
        """
        S_r = np.zeros(self.nRad_S)
        for i, surf in enumerate(self.radiatingSurface):
            surfPoints, group_indices = get_group_points(self.grid_show, surf)

            for j, f in enumerate(self.freq_array):
                self.Zab_total[j, i] = np.sum(self.p_total_mesh[j].coefficients[surfPoints]) / -self.surfaceVelocity[i][j]
                self.Zab_single[j, i] = np.sum(self.p_mesh[j, i].coefficients[surfPoints]) / -self.surfaceVelocity[i][j]
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

        micPosition = np.array(micPosition).T
        nMic = np.shape(micPosition)[1]

        pressure_mic_array = np.zeros([len(self.freq_array), nMic, self.nRad_S], dtype=complex)
        pressure_mic = np.zeros([len(self.freq_array), nMic], dtype=complex)
        omega = 2 * np.pi * self.freq_array
        k = - omega / self.c_0

        print("\n" + "Computing pressure at microphones")
        for i in tqdm(range(len(k))):  # looping through frequencies
            for rs in range(self.nRad_S):
                # print(i)
                # print(rs)
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

        pl.background_color = 'lightyellow'
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
                        n_colors=20, show_scalar_bar=True, cmap="turbo", scalar_bar_args=sargs, name='SPL',
                        lighting=False)

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

    def plot_admittanceRepartition(self, freq):
        grid = self.grid
        Yn = self.admittanceCoeff
        idx, _ = gtb.findInArray(self.freq_array, freq)
        return bempp.api.GridFunction(self.spaceP, coefficients=Yn[:, idx])


## ===========================
# %% preprocessing - acoustics
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
        vertex_location = centroids[slices[0], :]  # get [x, y, z] vertex position of radiator
        vertex_center = gtb.geometry.recenterZero(vertex_location)  # recenter geometry at [x=0, y=0, z=0]
    else:
        vertex_location = centroids[support_elements, :]  # get [x, y, z] vertex position of radiator
        vertex_center = gtb.geometry.recenterZero(vertex_location)  # recenter geometry at [x=0, y=0, z=0]

    ind = np.zeros(len(vertex_center))
    val = np.zeros([len(vertex_center), 3])
    ind_vert = []
    for i in range(len(vertex_center)):
        r = vertex_center[i, :]
        index, value = gtb.geometry.findClosestPoint(vibrometry_points, r[0], r[1], r[2])
        ind[i] = index  # Rc[index] that corresponds to the closest point to vertice_location_center
        val[i, :] = value

    Nfft = np.shape(vibrometric_data)[1]
    coefficients = np.zeros([Nfft, len(vertex_center)], dtype=complex)
    for i in range(len(vertex_center)):
        # numPoints = len(np.argwhere(ind == ind[i]))
        coefficients[:, i] = vibrometric_data[int(ind[i]), :]

    if sizeFactor > 1:
        coefficients = np.tile(coefficients, (1, sizeFactor))
        # print(coefficients.shape)
    return coefficients

def getSurfaceAdmittance(absorbingSurface, surfaceImpedance, freq, spaceP, c_0, rho_0):
    """
    Compute the total single layer coefficients linked to surfaces impedance.
    :param absorbingSurface:
    :param surfaceImpedance:
    :param freq:
    :param spaceP:
    :return:
    """
    Nfft       = len(freq)
    k          = - 2*np.pi*freq / c_0
    absSurf_in = np.array(absorbingSurface)
    surfImp_in = surfaceImpedance
    Nsurf      = len(absSurf_in)             # number of absorbing surfaces
    grid       = spaceP.grid
    gridLength = len(grid.domain_indices)    # number of indices in grid (domain of study)
    dofCount   = spaceP.grid_dof_count

    admittanceMatrix = np.ones([dofCount, Nfft], dtype=complex) * 0.0256 # default 5% damping
    if absSurf_in.shape[0] == 0 :
        None
    else:
        for surf in range(Nsurf):
            tmp_surf = absSurf_in[surf]  # current surface on which we apply admittance coefficients
            vertex, _ = get_group_points(grid, tmp_surf)
            for f in range(Nfft):
                try:
                    Yn = rho_0 * c_0 / surfImp_in[surf][f]
                except:
                    Yn = rho_0 * c_0 / surfImp_in[surf]
                admittanceMatrix[vertex, f] = np.ones(len(vertex)) * Yn
    return admittanceMatrix


## =============================
# %% useful functions - geometry
def get_group_points(grid, group_number):
    domain_indices = grid.domain_indices
    elements       = grid.elements
    element_edges  = grid.element_edges
    edges          = grid.edges
    # Find the indices where domain_indices match the group_number
    group_indices = np.where(domain_indices == group_number)  # indices of support elements

    # Get the corresponding columns of elements
    group_points = elements[:, group_indices]  # segments (and NOT points!) part of group_number
    # group_points   = element_edges[:, group_segments] # points part of group_number

    # Reshape to a 1D array and use unique to get unique values
    unique_points   = np.unique(group_points)

    return unique_points, group_indices

def getTriangleSurface(vertices, column_index):
    # Extract the three vertices of the triangle from the specified column
    triangle = vertices[:, column_index]

    # Unpack the vertices
    x1, y1, z1 = triangle[0]
    x2, y2, z2 = triangle[1]
    x3, y3, z3 = triangle[2]

    # Calculate the sides of the triangle
    side1 = np.linalg.norm(np.array([x2, y2, z2]) - np.array([x1, y1, z1]))
    side2 = np.linalg.norm(np.array([x3, y3, z3]) - np.array([x2, y2, z2]))
    side3 = np.linalg.norm(np.array([x1, y1, z1]) - np.array([x3, y3, z3]))

    # Calculate the semi-perimeter
    s = (side1 + side2 + side3) / 2

    # Calculate the area using Heron's formula
    area = np.sqrt(s * (s - side1) * (s - side2) * (s - side3))
    return area