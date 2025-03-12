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
from ..general import plot as gplot
directory_path = os.path.abspath(electroacPy.__file__)

pi = np.pi

class evaluations:
    """
    Tools for creating and managing various types of acoustic evaluations 
    using a Boundary Element Method (BEM) object.

    Available evaluation types:
    - Planar evaluations: Defined on a planar surface with specified dimensions.
    - Polar evaluations: Defined along a circular path in a specified plane.
    - Spherical evaluations: Defined on a spherical surface.
    - Pressure response evaluations: Recorded pressure responses at specified microphone positions.
    - Custom plotting grids through import
    
    """
    def __init__(self, bemObject):
        self.bemObject = bemObject
        self.frequency = bemObject.frequency
        self.setup = {}
        
        # ref to system
        self.referenceStudy = None
    
    
    def polarRadiation(self, evaluationName, minAngle: float, maxAngle: float,
                       step: float, on_axis: str, direction: str, 
                       radius: float = 5, offset: list = [0, 0, 0], **kwargs):
        
        if evaluationName not in self.setup:
            self.setup[evaluationName] = PolarRadiation(minAngle, maxAngle, 
                                                         step, on_axis, direction,
                                                         radius, offset)
        elif evaluationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(evaluationName)
            self.setup[evaluationName] = PolarRadiation(minAngle, maxAngle, 
                                                         step, on_axis, direction,
                                                         radius, offset)
        else:
             print("evaluation {} already exists. You can overwrite it using  \
                   the 'overwrite' flag".format(evaluationName))   
            
    
    def pressureField(self, evaluationName, Length: float, Width: float,
                      step: float, plane: str, offset: list = [0, 0, 0],
                      **kwargs):
        
        if evaluationName not in self.setup:
            self.setup[evaluationName] = PressureField(Length, Width, step, 
                                                        plane, offset)
        elif evaluationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(evaluationName)
            self.setup[evaluationName] = PressureField(Length, Width, step, 
                                                        plane, offset)
        else:
             print("evaluation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(evaluationName))

        
    def fieldPoint(self, evaluationName, microphonePosition, **kwargs):        
        if evaluationName not in self.setup:
            self.setup[evaluationName] = FieldPoint(microphonePosition, **kwargs)

        elif evaluationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(evaluationName)
            self.setup[evaluationName] = FieldPoint(microphonePosition, **kwargs)

        else:
             print("evaluation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(evaluationName))   


    def boundingBox(self, evaluationName, Lx, Ly, Lz,
                    step=1, offset=[0, 0, 0], **kwargs):
        if evaluationName not in self.setup:
            self.setup[evaluationName] = BoundingBox(Lx, Ly, Lz, step, offset)

        elif evaluationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(evaluationName)
            self.setup[evaluationName] = BoundingBox(Lx, Ly, Lz, step, offset)

        else:
             print("evaluation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(evaluationName))   
        
        
    def sphericalRadiation(self, evaluationName, nMic, 
                           radius, offset=[0, 0, 0], **kwargs):
        if evaluationName not in self.setup:
            self.setup[evaluationName] = SphericalRadiation(nMic, radius, offset)

        elif evaluationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(evaluationName)
            self.setup[evaluationName] = SphericalRadiation(nMic, radius, offset)

        else:
             print("evaluation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(evaluationName))  
                   
    def plottingGrid(self, evaluationName, path_to_grid, **kwargs):
        if evaluationName not in self.setup:
            self.setup[evaluationName] = PlottingGrid(path_to_grid)

        elif evaluationName in self.setup and "overwrite" in kwargs:
            self.setup.pop(evaluationName)
            self.setup[evaluationName] = PlottingGrid(path_to_grid)

        else:
             print("evaluation {} already exists. You can overwrite it using \
                   the 'overwrite' flag".format(evaluationName))  
    
    def solve(self, evaluation_name="all"):
        if evaluation_name == "all":
            obs_to_compute = [] #list(self.setup.keys())
            for key in self.setup:
                if self.setup[key].isComputed is False:
                    obs_to_compute.append(key)
                else:
                    pass
        elif isinstance(evaluation_name, str):
            obs_to_compute = [evaluation_name]
        elif isinstance(evaluation_name, list):
            obs_to_compute = evaluation_name
        else:
            raise ValueError("evaluation_name should be a key in setup.")
        
        # group microphones to compute in one step (faster for large numbers
        # of microphones)
        mic2compute = np.empty([0, 3])
        for obs in obs_to_compute:
            setup = self.setup[obs]
            xMic = setup.xMic
            nMic = setup.nMic
            if setup.isComputed is False:
                mic2compute = np.concatenate((mic2compute, xMic))
        
        if mic2compute.any() == np.True_:
            _, pMic = self.bemObject.getMicPressure(mic2compute, 
                                                    individualSpeakers=True)
        else:
            None
        
        # ungroup microphones and store within each setup
        current_index = 0
        for obs in obs_to_compute:
            setup = self.setup[obs]
            nMic = setup.nMic
            setup.pMic = pMic[:, current_index:current_index+nMic, :]
            setup.isComputed = True
            current_index += nMic

    def plot_system(self, backend="pyvista"):
        """
        Plot set of evaluations.

        Parameters
        ----------
        backend : str, optional
            Which software is used for display. Can be either "pyvista" or "gmsh".
            The default is "pyvista".

        Returns
        -------
        obj : None
            Plot.

        """
        if backend == "pyvista":
            obj = plot_system_pv(self)
        elif backend == "gmsh":
            obj = plot_system_gmsh(self)
        return obj
    
    def plot(self, evaluations=[], radiatingElement=[], processing=None,
             pf2grid=False, **kwargs):
        """
        Plot evaluations for given radiatingElement

        Parameters
        ----------
        evaluations : str or list of str, optional
            List of evaluations to plot. The default is "all".
        radiatingElement : int or list of int, optional
            List of radiating element to plot. The default is "all"
        processing : postProcessing object
            Post-processing to apply on given surfaces.
        pf2grid : bool
            If True, pressure-fields are "converted" into grids to plot in gmsh.
        **kwargs : key arguments for post-processing options
            for now there is not much options. Mostly for transformation of 
            plotting grids.
                - transformation : str
                    Either 'real', 'imag', 'spl' or 'SPL', 'phase'

        Returns
        -------
        None.

        """
        # Check input for which data to plot
        if bool(evaluations) is False:
            # plot all evaluations
            obs2plot = list(self.setup.keys())
        elif isinstance(evaluations, list):
            obs2plot = evaluations
        else:
            obs2plot = [evaluations]
        
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
        
        # Sort evaluations by type -> this could def be better, but it works sooo...
        polar, polarName = [], []
        field, pmicField, L, W, xMic_f = [], [], [], [], []
        point = []
        box = []
        sphere = []
        grid = []
        gridName = []
        for key in obs2plot:
            setup = self.setup[key]
            if setup.type == "polar":
                polar.append(setup)
                polarName.append(key)
            elif setup.type == "pressure_field" and pf2grid is False:
                field.append(setup)
                pmicField.append(setup.pMic)
                xMic_f.append(setup.xMic)
                L.append(setup.L)
                W.append(setup.W)
            elif setup.type == "box":
                box.append(setup)
            elif setup.type == "spherical":
                sphere.append(setup)
            elif setup.type == "field_point":
                point.append(setup)
            elif setup.type == "grid":
                grid.append(setup)
                gridName.append(key)
            elif setup.type == "pressure_field" and pf2grid is True:
                grid.append(setup)
                gridName.append(key)
        
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
            gplot.pressureField_3D(self.bemObject, 
                                   xMic_f, L, W, field2plot, element2plot)
        
        # plot evaluation points
        for i, obs in enumerate(point):
            pMic = obs.pMic
            point2plot = getPressure(pMic, self.bemObject.radiatingElement, 
                                    element2plot, elementCoeff)
            gplot.FRF(self.frequency, point2plot, legend=obs.legend)
            
        # plot spherical radiation
        for i, obs in enumerate(sphere):
            pMic = obs.pMic
            xMic = obs.xMic
            point2plot = getPressure(pMic, self.bemObject.radiatingElement,
                                     element2plot, elementCoeff)
            gplot.sphericalRadiation(xMic, point2plot, self.frequency)
            
        # plot bounding box
        for i, obs in enumerate(box):
            nx, ny, nz = obs.nx, obs.ny, obs.nz
            pMic = obs.pMic
            xMic = obs.xMic
            point2plot = getPressure(pMic, self.bemObject.radiatingElement,
                                     element2plot, elementCoeff)
            gplot.boundingBox(self.bemObject, nx, ny, nz, point2plot,
                              xMic, radiatingElement)
            
        if bool(grid) is True:
            point2plot = []

            for i, obs in enumerate(grid):
                pMic = obs.pMic
                point2plot.append(getPressure(pMic, self.bemObject.radiatingElement,
                                              element2plot, elementCoeff))
            
            if "transformation" in kwargs:
                transformation = kwargs["transformation"]
            else:
                transformation = "SPL"
            if "export_grid" in kwargs:
                export_grid = kwargs["export_grid"]
            else:
                export_grid = False
            gplot.bempp_grid(self.bemObject, grid, point2plot,
                             elementCoeff, element2plot, gridName, 
                             transformation, export_grid)
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
        
        from electroacPy.general.geometry import create_circular_array
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
        if isinstance(micPositions, np.ndarray):
            self.xMic = micPositions    
        else:
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
        
        from electroacPy.general.geometry import create_planar_array
        self.xMic, self.L, self.W = create_planar_array(Length, Width, 
                                                        step, plane, offset)        
        
        self.nMic = len(self.xMic)
            
        # help with dimensions - might need to change
        self.nx, self.ny, self.nz = 0, 0, 0
        L, W = len(self.L), len(self.W)
        
        # Mapping planes to dimensions
        plane_mapping = {
            ('x', 'y'): (L, W, 1),
            ('x', 'z'): (L, 1, W),
            ('y', 'x'): (W, L, 1),
            ('y', 'z'): (1, L, W),
            ('z', 'x'): (W, 1, L),
            ('z', 'y'): (1, W, L)
        }
        
        # Set dimensions based on the selected plane
        self.nx, self.ny, self.nz = plane_mapping.get((self.plane[0], 
                                                      self.plane[1]), 
                                                      (0, 0, 0))

        # edit evaluationParameters
        self.type = "pressure_field"
        self.isComputed = False
        self.pMic = None
        
        # data for gmsh plots
        from scipy.spatial import Delaunay
        if plane == "xy" or plane == "yx":
            dim_ = [0, 1]
        elif plane == "xz" or plane == "zx":
            dim_ = [0, 2]
        elif plane == "yz" or plane == "zy":
            dim_ = [1, 2]
    
        triangulation = Delaunay(self.xMic[:, dim_])
        self.elements = triangulation.simplices.T
        self.dof      = self.nMic
        self.vertices = self.xMic.T

class BoundingBox:
    def __init__(self, Lx, Ly, Lz, step=1, offset=[0, 0, 0]):
        self.Lx = Lx
        self.Ly = Ly 
        self.Lz = Lz
        self.step = step
        self.offset = offset
        
        from electroacPy.general.geometry import create_bounding_box
        self.xMic, dim = create_bounding_box(Lx, Ly, Lz, step, offset)
        self.nMic = len(self.xMic)
        self.nx = dim[0]
        self.ny = dim[1]
        self.nz = dim[2]

        # edit evaluationParameters
        self.type = "box"
        self.isComputed = False
        self.pMic = None
        
        
class SphericalRadiation:
    def __init__(self, nMic, radius, offset):
        self.radius = radius
        self.offset = offset
        
        from electroacPy.general.geometry import create_spherical_array
        self.xMic = create_spherical_array(nMic, radius, offset)
        self.nMic = len(self.xMic)

        # edit evaluationParameters
        self.type = "spherical"
        self.isComputed = False
        self.pMic = None
    

class PlottingGrid:
    def __init__(self, path_to_grid):
        from bempp.api import import_grid
        from ..general.geometry import check_mesh
        
        self.path_to_grid = check_mesh(path_to_grid)
        grid = import_grid(self.path_to_grid)
        self.dof      = grid.vertices.shape[1]
        self.vertices = grid.vertices
        self.elements = grid.elements
        self.edges    = grid.edges[:, grid.edge_on_boundary]
        # self.edges_on_boundary = grid.edge_on_boundary
        self.xMic     = grid.vertices.T
        self.nMic     = len(self.xMic)
        
        # edit evaluationParameters
        self.type = "grid"
        self.isComputed = False
        self.pMic = None
        
#%% Plotting help

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
    planes = pyvista.PolyData()
    
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
            
            plane = pyvista.Plane(center=point, direction=normal, 
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
            
            plane = pyvista.Plane(center=point, direction=normal, 
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
            
            plane = pyvista.Plane(center=point, direction=normal, 
                             i_size=plane_size[0]-x_offset,
                             j_size=plane_size[1]-y_offset)
        
        else:
            raise ValueError(f"Unknown boundary axis: {axis}. Must be 'x', 'y', or 'z'.")
        
        # Append the generated plane to the overall PolyData object
        planes += plane
    
    return planes
        
 
#%% PLOT SYSTEM
def plot_system_pv(eval_obj):
    """
    Plot current evaluation points in PyVista.
    
    Parameters
    ----------
    eval_obj: evaluation object
        Display set of evaluation.
    """
    
    if hasattr(eval_obj.bemObject, "xSource") and eval_obj.bemObject.meshPath is None:
        pl = pyvista.Plotter()
        pl.disable_anti_aliasing()
        xSource = eval_obj.bemObject.xSource
        point_cloud_source = pyvista.PolyData(xSource.astype(float))
        pl.add_mesh(point_cloud_source, color="black",
                    render_points_as_spheres=True, label="Sources", point_size=6)
        boundary_parameters = eval_obj.bemObject.boundary_conditions.parameters
    else:
        mesh = pyvista.read(eval_obj.bemObject.meshPath)
        radSurf = eval_obj.bemObject.radiatingElement
        boundary_parameters = eval_obj.bemObject.boundary_conditions
        
        # prepare mesh color, camera and bounds
        colors = []
        for scalar in mesh.active_scalars:
            if scalar in radSurf:
                colors.append([255, 0, 0])  # Red color for radiators
            else:
                colors.append([128, 128, 128])  # gray color for non-radiating
        mesh.cell_data['colors'] = colors
    
        # create plotter
        pl = pyvista.Plotter()
        pl.disable_anti_aliasing()
        if eval_obj.bemObject.domain == "exterior":
            pl.add_mesh(mesh, show_edges=True, cmap='summer', scalars='colors',
                        show_scalar_bar=False, rgb=True)
        elif eval_obj.bemObject.domain == "interior":
            pl.add_mesh_clip_box(mesh, show_edges=True, cmap='summer',  
                                 scalars='colors', show_scalar_bar=False, 
                                 rgb=True, rotation_enabled=False)
        light = pyvista.Light(light_type='headlight')
        pl.add_light(light)
    
    
    
    # add evaluation points
    colour = ['blue', 'red', 'green', 'orange', 
              'purple', 'brown', 'yellow', 'cyan'] * 4
    for i, key in enumerate(eval_obj.setup):
        setup = eval_obj.setup[key]
        xMic = setup.xMic
        point_cloud_tmp = pyvista.PolyData(xMic.astype(float))
        if setup.type == "box":
            pl.add_mesh(point_cloud_tmp.outline(), color=colour[i],
                        render_points_as_spheres=True, label=key, point_size=6)
        elif setup.type == "pressure_field":
            grid = pyvista.StructuredGrid()
            grid.points = xMic
            grid.dimensions = [setup.nx, setup.ny, setup.nz]
            pl.add_mesh(grid, show_edges=True, 
                        style='wireframe', color=colour[i], label=key)
        else:
            pl.add_mesh(point_cloud_tmp, color=colour[i],
                        render_points_as_spheres=True, label=key, point_size=6)

    # get bounds to add planes
    bounds = pl.bounds
    x_width = bounds[1] - bounds[0]  # Width along the x-axis
    y_width = bounds[3] - bounds[2]  # Width along the y-axis
    z_width = bounds[5] - bounds[4]  # Width along the z-axis
    
    if bool(eval_obj.setup) is False:   
        p_coeff = 2  # make infinite boundary a bit bigger if no eval is set 
    else:
        p_coeff = 1
    
    plane_size = (x_width*p_coeff, y_width*p_coeff, z_width*p_coeff)
    # add floor if infinite boundary conditions are present in bemObject
    bemObject = eval_obj.bemObject
    if bemObject.boundary_conditions is not None:
        inf_bc = []
        offset = []
        for key in boundary_parameters:
            if key in ["x", "X", "y", "Y", "z", "Z"]:
                inf_bc.append(key)
                offset.append(boundary_parameters[key]["offset"]) 
        
        if bool(inf_bc) is True:
            boundary_planes = create_boundary_planes(inf_bc, offset, plane_size)
            # Add boundary planes to plotter
            pl.add_mesh(boundary_planes, color="gray", opacity=0.5)
        else:
            pass

    if bool(eval_obj.setup) is True:   
        pl.add_legend(face='circle', bcolor=None)

    pl.add_axes(color='black')
    pl.background_color = 'white'

    _ = pl.show_grid(color='k')
    obj = pl.show()
    return obj


def plot_system_gmsh(eval_obj):
    """
    Plot current evaluation points in GMSH.
    
    Parameters
    ----------
    eval_obj: evaluation object
        Display set of evaluation.
    """
    import gmsh
    gmsh.initialize()
    
    if hasattr(eval_obj.bemObject, "xSource") and eval_obj.bemObject.meshPath is None:
        source_group = []
        for i in range(len(eval_obj.bemObject.xSource)):
            xSource = eval_obj.bemObject.xSource[i, :]
            point = gmsh.model.geo.add_point(xSource[0], xSource[1], xSource[2])
            source_group.append(point) 
            gmsh.model.geo.synchronize()
            gmsh.model.setColor([(0, point)], 0, 0, 0)
        physical_group_id = gmsh.model.addPhysicalGroup(0, source_group)  # '0' stands for point entities
        gmsh.model.setPhysicalName(0, physical_group_id, "sources")  # Optional: name the physical group
    elif hasattr(eval_obj.bemObject, "xSource") and eval_obj.bemObject.meshPath is not None:
        gmsh.open(eval_obj.bemObject.meshPath)
        source_group = []
        for i in range(len(eval_obj.bemObject.xSource)):
            xSource = eval_obj.bemObject.xSource[i, :]
            point = gmsh.model.geo.add_point(xSource[0], xSource[1], xSource[2])
            source_group.append(point) 
            gmsh.model.geo.synchronize()
            gmsh.model.setColor([(0, point)], 0, 0, 0)
        physical_group_id = gmsh.model.addPhysicalGroup(0, source_group)  # '0' stands for point entities
        gmsh.model.setPhysicalName(0, physical_group_id, "sources")  # Optional: name the physical group
    else:
        gmsh.open(eval_obj.bemObject.meshPath)
  

    # get max nodes in system
    all_nodes = gmsh.model.mesh.getNodes()
    maxNode =  []
    for i in all_nodes:
        try:
            maxNode.append(np.max(i))
        except:
            None
    maxNODE = int(np.max(maxNode)) + 1
    
    colors = [
        [0, 0, 255],    # Blue
        [255, 165, 0],  # Orange
        [0, 128, 0],    # Green
        [255, 0, 0],    # Red
        [128, 0, 128],  # Purple
        [165, 42, 42],  # Brown
        [255, 192, 203],# Pink
        [128, 128, 128],# Grey
        [128, 128, 0],  # Olive
        [0, 255, 255],  # Cyan
    ]

    # Function to repeat the colors to match the number of observations
    def get_colors(num_observations):
        return [colors[i % len(colors)] for i in range(num_observations)]

    nEval = len(eval_obj.setup)
    colorsEval = get_colors(nEval) 

    for s, setup in enumerate(eval_obj.setup):
        cl = colorsEval[s]
        csetup = eval_obj.setup[setup]
        obs_group = []
        nMic = csetup.nMic
        xMic = csetup.xMic
        if csetup.type not in ["grid", "pressure_field"]:
            for i in range(nMic):
                point = gmsh.model.geo.add_point(xMic[i, 0], xMic[i, 1], xMic[i, 2])
                obs_group.append(point) 
                gmsh.model.geo.synchronize()
                gmsh.model.setColor([(0, point)], cl[0], cl[1], cl[2])     # Red for point p1
            physical_group_id = gmsh.model.addPhysicalGroup(0, obs_group)  # '0' stands for point entities
            gmsh.model.setPhysicalName(0, physical_group_id, setup)        # Optional: name the physical group
        
        elif csetup.type in ["grid", "pressure_field"]:
            # node
            nodeTags_n = []
            coord = []
            
            try: # in case plottingGrids are set before other evaluations
                maxDim0 = int(np.max(gmsh.model.getEntities(0)) + 1)
            except:
                maxDim0 = 1    
            maxDim2 = int(np.max(gmsh.model.getEntities(2)) + 1)
            
        
            # elements
            elements              = csetup.elements
            nElements             = elements.shape[1]
            nodeTags_elements     = []
            elementTags_elements  = []
            
            for i in range(nMic):
                nodeTags_n.append(i+maxNODE)
                coord.append(xMic[i, 0])
                coord.append(xMic[i, 1])
                coord.append(xMic[i, 2])
        
            # elements of dim 2
            for k in range(nElements):
                elementTags_elements.append(k+1)
                nodeTags_elements.append(elements[0, k]+maxNODE) 
                nodeTags_elements.append(elements[1, k]+maxNODE) 
                nodeTags_elements.append(elements[2, k]+maxNODE)
                

            gmsh.model.addDiscreteEntity(0, maxDim0)
            gmsh.model.addDiscreteEntity(2, maxDim2)
            gmsh.model.mesh.addNodes(0, maxDim0, nodeTags_n, coord)
            gmsh.model.mesh.addElements(2, maxDim2, [2], 
                                                    [elementTags_elements], [nodeTags_elements])
            
            physical_group_id = gmsh.model.addPhysicalGroup(2, [maxDim2])  # '0' stands for point entities
            gmsh.model.setPhysicalName(2, physical_group_id, setup) 
            
            maxNODE += i
    
    gmsh.option.setNumber("Geometry.PointSize", 3)
    gmsh.option.setNumber("Geometry.PointType", 1)
    gmsh.model.geo.synchronize()

    gmsh.fltk.run()

    gmsh.finalize()
    return None



def create_gmsh_mesh(xMic, elements, mesh_filename="mesh.msh"):
    import gmsh
    gmsh.initialize()
    gmsh.model.add("mic_mesh")

    # Add nodes
    nMic = xMic.shape[0]
    for i in range(nMic):
        gmsh.model.mesh.addNode(i + 1, xMic[i, 0], xMic[i, 1], xMic[i, 2])

    # Add triangular elements
    nElements = elements.shape[1]
    elementType = gmsh.model.mesh.getElementType("triangle", 1)  # 2D triangle element
    gmsh.model.mesh.addElements(2, [elementType], [[i + 1 for i in range(nElements)]], elements + 1)