"""
Functions to manipulate geometrical data (point cloud, etc.)
"""

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.spatial import Delaunay


def check_mesh(mesh_path):  
    """
    Check input format of a given mesh. Only .msh (v2 ASCII) and .med files 
    are currently supported.
    
    """
    if mesh_path[-4:] == ".msh":
        meshFile = open(mesh_path)
        lines = meshFile.readlines()
        if lines[1][0] != '2':
            raise TypeError(
                "Mesh file is not in version 2. Errors will appear when mirroring mesh along boundaries.")
        meshFile.close()
        mesh_path_update = mesh_path
        
    elif mesh_path[-4:] == ".med": # conversion from med to msh to keep groups
        import gmsh
        print("\n")
        print("Conversion from *.med to *.msh... \n")
        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.open(mesh_path)
        gmsh.write(mesh_path[:-4] + ".msh")
        gmsh.finalize()
        mesh_path_update = mesh_path[:-4] + ".msh"
        
    else: # conversion from med to msh to keep groups
        mesh_path_update = mesh_path
        raise Exception(
            "Not compatible file format. Try *.med or *.msh.")
    return mesh_path_update

def extract_bem_surface(grid):
    """
    Return index of nodes corresponding grouped surface.
    Extract from a bempp grid
    :param grid: 
    :return: 
    """
    vertices = grid.vertices
    elements = grid.elements
    domain_indices = grid.domain_indices
    domains = np.unique(domain_indices)
    groups = {}
    nodes = {}
    for j in domains:
        # initialise groups
        groups[str(j)] = []
        nodes[str(j)] = None

    for i in range(len(domain_indices)):
        for j in domains:
            if j == domain_indices[i]:
                groups[str(j)] += list(elements[:, i])

    for j in domains:
        nodes[str(j)] = np.unique(groups[str(j)])
    
    return nodes


def findClosestPoint(reference_array, x, y, z):
    """
    Finds the index and coordinates of the closest point in the reference array
    to the specified coordinates (x, y, z).

    Parameters
    ----------
    reference_array : numpy array, shape (Npoints, 3)
        Array of reference points with shape [Npoints, 3].

    x, y, z : float
        Coordinates (x, y, z) to find the closest point to (in the reference array).

    Returns
    -------
    closest_index : int
        Index of the closest point in the reference array.

    closest_coordinates : numpy array, shape (3,)
        Coordinates of the closest point in the reference array.
    """
    # print(reference_array.shape)
    # print(reference_array[20, 0])

    # distances = np.linalg.norm(reference_array - np.array([x, y, z]), axis=1)
    distances = np.zeros(len(reference_array))
    for p in range(len(reference_array)):
        distances[p] = np.sqrt((reference_array[p, 0] - x)**2 + (reference_array[p, 1] - y)**2 +
                            (reference_array[p, 2] - z)**2)

    closest_index = np.argmin(distances)
    closest_coordinates = reference_array[closest_index]
    return closest_index, closest_coordinates


def points_within_radius(R_ref, R_meas, radius, v_meas, Nfft):
    """
    Find points in R_meas within a specified radius of points in R_ref.

    Parameters:
        R_ref: Array of reference points of shape (Npoints_ref, 3).
        R_meas: Array of measured points of shape (Npoints_meas, 3).
        radius: Radius within which to search for points.

    Returns:
        List of indices of points in R_meas within the specified radius of each point in R_ref.
    """
    
    indices_within_radius = []
    coefficients = np.zeros([len(R_ref), Nfft], dtype=complex)

    for ind_point, ref_point in enumerate(R_ref):
        inc = 0
        # Calculate distance between ref_point and all points in R_meas
        distances = np.zeros(len(R_meas))
        for i in range(len(R_meas)):
            distances[i] = np.sqrt((R_meas[i, 0] - R_ref[ind_point, 0])**2
                                   + (R_meas[i, 1] - R_ref[ind_point, 1])**2
                                   + (R_meas[i, 2] - R_ref[ind_point, 2])**2)

        # Find indices of points within radius - calculate corresponding weight
        within_radius_indices = np.array([])
        while len(within_radius_indices) == 0:
            within_radius_indices = np.where(distances <= radius + inc)[0]
            inc += 5e-4
        weight = 1 - distances[within_radius_indices] / radius
        weighted_average = np.sum(weight)

        # Add indices to the list
        indices_within_radius.append(within_radius_indices)

        # store coeff
        coefficients[ind_point, :] = np.sum(v_meas[within_radius_indices, :] * weight[:, np.newaxis],
                                            axis=0) / weighted_average

    return indices_within_radius, coefficients.T


def recenterZero(X):
    """
    Recenter a point cloud at the origin (0, 0, 0) in Cartesian coordinates.

    Parameters
    ----------
    X : numpy array
        Array of shape (nPoints, 3) representing the point cloud.

    Returns
    -------
    Xcenter : numpy array
        Point cloud centered at the origin.

    Notes
    -----
    This function takes a point cloud represented by a numpy array 'X' of shape (nPoints, 3),
    where each row corresponds to a 3D point with Cartesian coordinates (x, y, z). The function
    calculates the mean values of the x, y, and z coordinates and recenters the entire point cloud
    such that its center is at the origin (0, 0, 0). The resulting 'Xcenter' array contains the
    recentered point cloud.
    """
    mean_x = np.mean(X[:, 0])
    mean_y = np.mean(X[:, 1])
    mean_z = np.mean(X[:, 2])
    Xcenter = np.zeros([len(X), 3])
    Xcenter[:, 0] = X[:, 0] - mean_x
    Xcenter[:, 1] = X[:, 1] - mean_y
    Xcenter[:, 2] = X[:, 2] - mean_z
    return Xcenter

def rotatePointCloud(point_cloud, rotation_angles_deg):
    """
    Rotate a point cloud around x, y, and z axes.

    Parameters
    ----------
    point_cloud : numpy array, shape (nPoints, 3)
        Point cloud coordinates in the shape (x, y, z).

    rotation_angles_deg : list, shape (3,)
        Rotation angles around x, y, and z axes in degrees.

    Returns
    -------
    rotated_point_cloud : numpy array, shape (nPoints, 3)
        Rotated point cloud coordinates after applying the rotations.

    Notes
    -----
    This function rotates a given point cloud by specified angles around the x, y, and z axes.
    The input 'point_cloud' should be a numpy array with shape (nPoints, 3), representing the
    Cartesian coordinates (x, y, z) of each point in the cloud. The 'rotation_angles_deg' list
    should contain three rotation angles in degrees corresponding to the rotations around the
    x, y, and z axes, respectively. The rotations are applied sequentially in the order: x, y, z.
    The function returns the rotated point cloud coordinates as a numpy array with the same shape
    as the input 'point_cloud'.
    """

    # Convert rotation angles from degrees to radians
    rotation_angles_deg = np.array(rotation_angles_deg)
    rotation_angles_rad = np.radians(rotation_angles_deg)

    # Create rotation object for each axis
    rotation_x = Rotation.from_euler('x', rotation_angles_rad[0], degrees=False)
    rotation_y = Rotation.from_euler('y', rotation_angles_rad[1], degrees=False)
    rotation_z = Rotation.from_euler('z', rotation_angles_rad[2], degrees=False)

    # Apply rotations to the point cloud
    rotated_point_cloud = rotation_x.apply(point_cloud)
    rotated_point_cloud = rotation_y.apply(rotated_point_cloud)
    rotated_point_cloud = rotation_z.apply(rotated_point_cloud)

    return rotated_point_cloud


def compute_circle_point_cloud_surface_area(point_cloud):
    # Ensure that the point cloud is a numpy array of shape [nPoints, 3]
    if point_cloud.shape[1] != 3:
        raise ValueError("Point cloud should have shape [nPoints, 3]")

    # Perform Delaunay triangulation
    tri = Delaunay(point_cloud)
    maxR = np.max(tri.max_bound)
    surface_area = np.pi * maxR**2
    return surface_area


## MICROPHONE ARRAY CREATION
def create_circular_array(theta, on_axis, rotation, radius, offset):
    # Initialize the points array
    points = np.zeros((len(theta), 3))

    # Convert theta to radians
    theta_rad = np.deg2rad(theta)
    plane_indices = ()

    if on_axis in ["+x", "x", "-x"]:
        plane_indices += (0,)
    elif on_axis in ["+y", "y", "-y"]:
        plane_indices += (1,)
    elif on_axis in ["+z", "z", "-z"]:
        plane_indices += (2,)

    if rotation in ["+x", "x", "-x"]:
        plane_indices += (0,)
    elif rotation in ["+y", "y", "-y"]:
        plane_indices += (1,)
    elif rotation in ["+z", "z", "-z"]:
        plane_indices += (2,)

    # Determine the sign based on `on_axis`
    sign = 1 if "-" not in on_axis else -1
    sign_rotation = 1 if "-" not in rotation else -1

    if plane_indices[0] == plane_indices[1]:
        raise Exception("Direction of rotation cannot be similar to on-axis value.")

    # Generate the circular coordinates
    points[:, plane_indices[0]] = radius * np.cos(theta_rad) * sign
    points[:, plane_indices[1]] = radius * np.sin(theta_rad) * sign_rotation

    # Apply the offsets
    for i in range(3):
        points[:, i] += offset[i]

    return points


def create_planar_array(length, width, micSpacing, plane, offset=[0, 0, 0], 
                      vert=False, mode=False):
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
    nMic_L = int(length / micSpacing)
    nMic_W = int(width / micSpacing)
    L = np.linspace(0, length, nMic_L) #np.arange(0, length+micSpacing, micSpacing)
    W = np.linspace(0, width, nMic_W)  #np.arange(0, width+micSpacing, micSpacing)
    nMic = nMic_L * nMic_W  #len(L)*len(W)
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
        xmic_n = filter_points(xmic, vert, mode=mode)
        L_n = L[np.isin(L, xmic_n[:, dim1])]
        W_n = W[np.isin(W, xmic_n[:, dim2])]
        out = (xmic_n, L_n, W_n)
    return out


def filter_points(points, boundary, mode='inside'):
    """
    Filter points based on their location relative to the boundary.

    Parameters:
    - microphones: NumPy array of shape (N, 3), representing points coordinates.
    - boundary: NumPy array of shape (M, 3), representing boundary coordinates.
    - mode: 'inside' or 'outside', specifying whether to keep microphones inside or outside the boundary.

    Returns:
    - filtered_points: NumPy array of shape (K, 3), where K <= N, containing the filtered microphone coordinates.
    """

    if mode not in ('inside', 'outside'):
        raise ValueError("Mode must be 'inside' or 'outside'.")

    if mode == 'inside':
        # Keep microphones inside the boundary
        condition = np.all((points[:, np.newaxis] >= boundary.min(axis=0)) &
                           (points[:, np.newaxis] <= boundary.max(axis=0)), axis=2)
    else:
        # Keep microphones outside the boundary
        condition = np.any((points[:, np.newaxis] < boundary.min(axis=0)) |
                           (points[:, np.newaxis] > boundary.max(axis=0)), axis=2)

    filtered_points = points[condition.all(axis=1)]
    return filtered_points

def create_spherical_array(nMic, sphereRadius=1.8, offset=[0, 0, 0]):
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


def create_bounding_box(Lx, Ly, Lz, step=1, offset=[0, 0, 0]):
    x_offset, y_offset, z_offset = offset
    x_range = np.arange(x_offset, x_offset + Lx + step, step)
    y_range = np.arange(y_offset, y_offset + Ly + step, step)
    z_range = np.arange(z_offset, z_offset + Lz + step, step)
    
    nx = len(x_range)
    ny = len(y_range)
    nz = len(z_range)
    dim = [nx, ny, nz]
    
    x_points, y_points, z_points = np.meshgrid(x_range, y_range, z_range)

    points = np.vstack([x_points.flatten(), y_points.flatten(), z_points.flatten()]).T

    return points, dim