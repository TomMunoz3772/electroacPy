#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:29:39 2024

@author: tom
"""

import numpy as np
import pyvista
import os
import electroacPy
import generalToolbox.plot as gplot
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
        self.frequency = bemObject.frequency
        self.setup = {}
        
        # ref to system
        self.referenceStudy = None
    
    
    def polarRadiation(self, observationName, minAngle: float, maxAngle: float,
                       step: float, on_axis: str, direction: str, 
                       radius: float = 5, offset: list = [0, 0, 0], **kwargs):
        
        if observationName not in self.setup:
            self.setup[observationName] = PolarRadiation(minAngle, maxAngle, 
                                                         step, on_axis, direction,
                                                         radius, offset)
        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = PolarRadiation(minAngle, maxAngle, 
                                                         step, on_axis, direction,
                                                         radius, offset)
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
            self.setup[observationName] = FieldPoint(microphonePosition, **kwargs)

        elif observationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(observationName)
            self.setup[observationName] = FieldPoint(microphonePosition, **kwargs)

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
    
    def solve(self, observation_name="all"):
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
            nMic = setup.nMic
            setup.pMic = pMic[:, current_index:current_index+nMic, :]
            setup.isComputed = True
            current_index += nMic

    def plot_system(self):
        mesh = pyvista.read(self.bemObject.meshPath)
        radSurf = self.bemObject.radiatingElement
        
        # prepare mesh color, camera and bounds
        colors = []
        for scalar in mesh.active_scalars:
            if scalar in radSurf:
                colors.append([255, 0, 0])  # Red color for radiators
            else:
                colors.append([128, 128, 128])  # gray color for non-radiating
        mesh.cell_data['colors'] = colors
        center = mesh.center

        # create plotter
        pl = pyvista.Plotter()
        pl.add_mesh(mesh, show_edges=True, cmap='summer',  scalars='colors',
                    show_scalar_bar=False, rgb=True)
        light = pyvista.Light(light_type='headlight')
        pl.add_light(light)
        pl.camera.focal_point = center

        
        # add evaluation points
        colour = ['blue', 'red', 'green', 'orange', 
                  'purple', 'brown', 'yellow', 'cyan'] * 4
        for i, key in enumerate(self.setup):
            setup = self.setup[key]
            xMic = setup.xMic
            point_cloud_tmp = pyvista.PolyData(xMic.astype(float))
            if setup.type == "box":
                pl.add_mesh(point_cloud_tmp.outline(), color=colour[i],
                            render_points_as_spheres=True, label=key, point_size=6)
            else:
                pl.add_mesh(point_cloud_tmp, color=colour[i],
                            render_points_as_spheres=True, label=key, point_size=6)

        # get bounds to add planes
        bounds = pl.bounds
        x_width = bounds[1] - bounds[0]  # Width along the x-axis
        y_width = bounds[3] - bounds[2]  # Width along the y-axis
        z_width = bounds[5] - bounds[4]  # Width along the z-axis
        plane_size = (x_width, y_width, z_width)
        # add floor if infinite boundary conditions are present in bemObject
        bemObject = self.bemObject
        if bemObject.boundary_conditions is not None:
            inf_bc = []
            offset = []
            for key in bemObject.boundary_conditions:
                if key in ["x", "X", "y", "Y", "z", "Z"]:
                    inf_bc.append(key)
                    offset.append(bemObject.boundary_conditions[key]["offset"]) 
            
            if bool(inf_bc) is True:
                boundary_planes = create_boundary_planes(inf_bc, offset, plane_size)
                # Add boundary planes to plotter
                pl.add_mesh(boundary_planes, color="gray", opacity=0.5)
            else:
                pass

        pl.add_axes(color='black')
        pl.add_legend(face='circle', bcolor=None)
        pl.background_color = 'white'

        _ = pl.show_grid(color='k')
        obj = pl.show()
        return obj
    
    def plot(self, evaluations=[], radiatingElement=[], processing=None):
        """
        Plot evaluations for given radiatingElement

        Parameters
        ----------
        evaluations : str or list of str, optional
            List of evaluations to plot. The default is "all".
        radiatingElement : int or list of int, optional
            List of radiating element to plot. The default is "all"
        processing : postProcessing object
            post-processing to apply on given surfaces.

        Returns
        -------
        None.

        """
        # Check input for which data to plot
        if bool(evaluations) is False:
            # plot all observations
            obs2plot = list(self.setup.keys())
        else:
            obs2plot = evaluations
        
        if bool(radiatingElement) is False:
            # plot all elements
            element2plot = self.bemObject.radiatingElement
        else:
            element2plot = radiatingElement
        
        if processing is not None:
            elementCoeff = np.ones([len(self.frequency), len(element2plot)], 
                                   dtype=complex)
            pp = processing
            for name in pp.TF:
                for idx, element in enumerate(element2plot):
                    if element in pp.TF[name]["radiatingElement"]:
                        elementCoeff[:, idx] *= pp.TF[name]["H"]
                        
        else:
            elementCoeff = np.ones([len(self.frequency), len(element2plot)], 
                                   dtype=complex)
        
        # Sort observations by type
        polar, polarName = [], []
        field, pmicField, L, W, xMic = [], [], [], [], []
        point, pmicPoint = [], []
        box = []
        sphere = []
        for key in obs2plot:
            setup = self.setup[key]
            if setup.type == "polar":
                polar.append(setup)
                polarName.append(key)
            elif setup.type == "pressure_field":
                field.append(setup)
                pmicField.append(setup.pMic)
                xMic.append(setup.xMic)
                L.append(setup.L)
                W.append(setup.W)
            elif setup.type == "box":
                box.append(setup)
            elif setup.type == "sphere":
                sphere.append(setup)
            elif setup.type == "field_point":
                point.append(setup)
                pmicPoint.append(setup.pMic)
        
        # plot polar
        for i, obs in enumerate(polar):
            pmic2plot = getPressure(obs.pMic, self.bemObject.radiatingElement, 
                                    element2plot, elementCoeff)
            gplot.directivityViewer(obs.theta, self.frequency, pmic2plot,
                                    title=polarName[i])
        
        # plot field
        if bool(field) is True:
            field2plot = getPressure(pmicField, self.bemObject.radiatingElement, 
                                    element2plot, elementCoeff)
            gplot.pressureField_3D(self.bemObject, xMic, L, W, field2plot, element2plot)
        
        # plot evaluation points
        for i, obs in enumerate(point):
            point2plot = getPressure(pmicPoint[i], self.bemObject.radiatingElement, 
                                    element2plot, elementCoeff)
            gplot.FRF(self.frequency, point2plot, legend=obs.legend)
        return None
    


#%% helper function
def getPressure(pmic, radiatingElement, element2plot, coefficients):
    """
    get the pressure to plot given element2plot and their corresponding 
    coefficients

    Parameters
    ----------
    pmic : numpy array
        pressure returned by one evaluation.
    radiatingElement : list of int
        radiatingElement of the BEM object. Used to associate pMic to the 
        correct coefficients
    element2plot : list of int
        specific radiating element to plot.
    coefficients : numpy array
        coefficients to apply to corresponding radiating surfaces.

    Returns
    -------
    pmic_out : TYPE
        DESCRIPTION.

    """
    if isinstance(pmic, list):
        pmic_out = []
        for p in range(len(pmic)):
            nFreq, nMic, nRad = pmic[p].shape
            pmic_tmp = np.zeros([nFreq, nMic], dtype=complex)
            for i, e in enumerate(element2plot):
                for mic in range(nMic):
                    radE = np.argwhere(e == np.array(radiatingElement))[0][0]
                    pmic_tmp[:, mic] += pmic[p][:, mic, radE] * coefficients[:, i]
            pmic_out.append(pmic_tmp)
    else:      
        nFreq, nMic, nRad = pmic.shape
        pmic_out = np.zeros([nFreq, nMic], dtype=complex)
        for i, e in enumerate(element2plot):
            for mic in range(nMic):
                radE = np.argwhere(e == np.array(radiatingElement))[0][0]
                pmic_out[:, mic] += pmic[:, mic, radE] * coefficients[:, i]
    return pmic_out
            

#%% Evaluation classes
class PolarRadiation:
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
        self.isComputed = False
        self.pMic = None
        
        

class FieldPoint:
    def __init__(self, micPositions, **kwargs): 
        self.xMic = np.array(micPositions)
        self.nMic = len(self.xMic)
        self.legend = None
        
        if "legend" in kwargs:
            self.legend = kwargs["legend"]
        
        # edit evaluationParameters
        self.type = "field_point"
        self.isComputed = False
        self.pMic = None
        
        
class PressureField:
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
        self.xMic, self.L, self.W = create_planar_array(Length, Width, 
                                                        step, plane, offset)        
        self.nMic = len(self.xMic)

        # edit evaluationParameters
        self.type = "pressure_field"
        self.isComputed = False
        self.pMic = None
        

class BoundingBox:
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
        self.type = "box"
        self.isComputed = False
        self.pMic = None
        
        
        
class SphericalRadiation:
    def __init__(self, nMic, radius, offset):
        self.radius = radius
        self.offset = offset
        
        from generalToolbox.geometry import create_spherical_array
        self.xMic = create_spherical_array(nMic, radius, offset)
        self.nMic = len(self.xMic)

        # edit evaluationParameters
        self.type = "spherical"
        self.isComputed = False
        self.pMic = None
    
#%% Plotting help
import pyvista as pv
import numpy as np

# def create_boundary_planes(boundary, offset, plane_size=1e6):
#     """
#     Create PyVista PolyData for boundaries in the plotter based on the provided normals and offsets.
    
#     Parameters:
#         boundary (list): List of boundary normals (e.g., ["x", "y", "z"]).
#         offset (list): List of offsets corresponding to each boundary normal.
#         plane_size (float): Size of the plane to be created, extending equally in both directions.
    
#     Returns:
#         pyvista.PolyData: PolyData object representing the boundary planes.
#     """
    
#     if len(boundary) != len(offset):
#         raise ValueError("The boundary and offset lists must have the same length.")
    
#     # Initialize an empty PolyData object to collect all the planes
#     planes = pv.PolyData()
    
#     # Loop through each boundary and create a corresponding plane
#     for i, axis in enumerate(boundary):
#         # Create the points for the plane depending on the axis
#         if axis == "x":
#             # Plane orthogonal to the x-axis (yz-plane)
#             point = [offset[i], 0, 0]  # Plane is located at x = offset[i]
#             normal = [1, 0, 0]  # Normal vector in x-direction
#             plane = pv.Plane(center=point, direction=normal, i_size=plane_size, j_size=plane_size)
        
#         elif axis == "y":
#             # Plane orthogonal to the y-axis (xz-plane)
#             point = [0, offset[i], 0]  # Plane is located at y = offset[i]
#             normal = [0, 1, 0]  # Normal vector in y-direction
#             plane = pv.Plane(center=point, direction=normal, i_size=plane_size, j_size=plane_size)
        
#         elif axis == "z":
#             # Plane orthogonal to the z-axis (xy-plane)
#             point = [0, 0, offset[i]]  # Plane is located at z = offset[i]
#             normal = [0, 0, 1]  # Normal vector in z-direction
#             plane = pv.Plane(center=point, direction=normal, i_size=plane_size, j_size=plane_size)
        
#         else:
#             raise ValueError(f"Unknown boundary axis: {axis}. Must be 'x', 'y', or 'z'.")
        
#         # Append the generated plane to the overall PolyData object
#         planes += plane
    
#     return planes

def create_boundary_planes(boundary, offset, plane_size=10):
    """
    Create PyVista PolyData for boundaries with adjusted positioning to form corners.

    Parameters:
        boundary (list): List of boundary normals (e.g., ["x", "y", "z"]).
        offset (list): List of offsets corresponding to each boundary normal.
        plane_size (float): Maximum size of the plane to be created (e.g., 10m).

    Returns:
        pyvista.PolyData: PolyData object representing the boundary planes.
    """
    
    if len(boundary) != len(offset):
        raise ValueError("The boundary and offset lists must have the same length.")
    
    # Initialize an empty PolyData object to collect all the planes
    planes = pv.PolyData()
    
    # Define default plane extent (10m by 10m)    

    # Loop through each boundary and create a corresponding plane with adjusted size
    for i, axis in enumerate(boundary):
        # Create the points for the plane depending on the axis
        if axis in ["x", "X"]:
            # Plane orthogonal to the x-axis (yz-plane)
            x_center = offset[i]  # Set x at the offset
            
            if "y" in boundary:
                y_offset = offset[boundary.index("y")]
                y_center = (plane_size[1]-y_offset)/2 + y_offset
            elif "Y" in boundary:
                y_offset = offset[boundary.index("Y")]
                y_center = (plane_size[1]-y_offset)/2 + y_offset
            else:
                y_offset=0
                y_center=0
            
            if "z" in boundary:
                z_offset = offset[boundary.index("z")]
                z_center = (plane_size[2]-z_offset)/2 + z_offset
            elif "Z" in boundary:
                z_offset = offset[boundary.index("Z")]
                z_center = (plane_size[2]-z_offset)/2 + z_offset
            else:
                z_offset=0
                z_center=0

            point = [x_center, y_center, z_center]  # Plane's center
            normal = [1, 0, 0]  # Normal vector in x-direction
            
            plane = pv.Plane(center=point, direction=normal, 
                             i_size=plane_size[2]-z_offset, 
                             j_size=plane_size[1]-y_offset)
        
        elif axis in ["y", "Y"]:
            # Plane orthogonal to the y-axis (xz-plane)
            y_center = offset[i]  # Set y at the offset
            
            if "x" in boundary:
                x_offset = offset[boundary.index("x")]
                x_center = (plane_size[0]-x_offset)/2 + x_offset
            elif "X" in boundary:
                x_offset = offset[boundary.index("X")]
                x_center = (plane_size[0]-x_offset)/2 + x_offset
            else:
                x_offset=0
                x_center=0
            
            if "z" in boundary:
                z_offset = offset[boundary.index("z")]
                z_center = (plane_size[2]-z_offset)/2 + z_offset
            elif "Z" in boundary:
                z_offset = offset[boundary.index("Z")]
                z_center = (plane_size[2]-z_offset)/2 + z_offset
            else:
                z_offset=0
                z_center=0
            
            point = [x_center, y_center, z_center]  # Plane's center
            normal = [0, 1, 0]  # Normal vector in x-direction
            
            plane = pv.Plane(center=point, direction=normal, 
                             i_size=plane_size[2]-z_offset, 
                             j_size=plane_size[0]-x_offset)
        
        elif axis in ["z", "Z"]:
            # Plane orthogonal to the z-axis (xy-plane)
            z_center = offset[i]  # Set z at the offset
            
            if "x" in boundary:
                x_offset = offset[boundary.index("x")]
                x_center = (plane_size[0]-x_offset)/2 + x_offset
            elif "X" in boundary:
                x_offset = offset[boundary.index("X")]
                x_center = (plane_size[0]-x_offset)/2 + x_offset
            else:
                x_offset=0
                x_center=0
            
            if "y" in boundary:
                y_offset = offset[boundary.index("y")]
                y_center = (plane_size[1]-y_offset)/2 + y_offset
            elif "Y" in boundary:
                y_offset = offset[boundary.index("Y")]
                y_center = (plane_size[1]-y_offset)/2 + y_offset
            else:
                y_offset=0
                y_center=0
            
            point = [x_center, y_center, z_center]  # Plane's center
            normal = [0, 0, 1]  # Normal vector in x-direction
            
            plane = pv.Plane(center=point, direction=normal, 
                             i_size=plane_size[0]-x_offset,
                             j_size=plane_size[1]-y_offset)
        
        else:
            raise ValueError(f"Unknown boundary axis: {axis}. Must be 'x', 'y', or 'z'.")
        
        # Append the generated plane to the overall PolyData object
        planes += plane
    
    return planes
        
        