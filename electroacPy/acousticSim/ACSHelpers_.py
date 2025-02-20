"""
Helper functions for BEM and point-source methods

"""
import numpy as np
import bempp.api


#%% Green-function and monopole 
def incidenceCoeff(xSource, Grid, domain, dofs="centroids"):
    """
    Get cos(angle) of source incidence to grid normals. In combination to -1jk, will 
    give the normal derivative (not exactly sure if this is the right term).

    Parameters
    ----------
    xSource : ndarray
        Array of source position (Nsource, 3).
    Grid : bempp-cl grid object
        System grid.
    domain : str
        Domain of simulation. Either "interior" or "exterior".
    dofs : str, optional
        Which DOFs to select. The default is "centroids".

    Returns
    -------
    n0 : TYPE
        DESCRIPTION.

    """
    if domain == "interior":
        domain_operator = -1
    else:
        domain_operator = 1
    Nvert = Grid.vertices.shape[1]
    Norm  = domain_operator * Grid.normals.T
    Vert  = Grid.vertices
    Cent  = Grid.centroids.T
    Ncent = Grid.centroids.shape[0]
    Ns    = len(xSource)

    if dofs == "centroids":
        DOF = Cent
        n0 = np.zeros(Ncent)
        dir_vect = np.zeros([3, Ncent, Ns])
        nDOF = Ncent
    
    elif dofs == "vertices":
        DOF = Vert
        n0 = np.zeros(Nvert)
        dir_vect = np.zeros([3, Nvert, Ns])
        nDOF = Nvert
    
    for s in range(Ns):
        dir_vect[0, :, s] = DOF[0, :] - xSource[s, 0]
        dir_vect[1, :, s] = DOF[1, :] - xSource[s, 1]
        dir_vect[2, :, s] = DOF[2, :] - xSource[s, 2]
    
    for i in range(3):
        for s in range(Ns):
            norm_tmp = (dir_vect[0, i, s]**2 + 
                        dir_vect[1, i, s]**2 + dir_vect[2, i, s]**2)**0.5
            dir_vect[:, i, s] /= norm_tmp #np.max(np.abs(dir_vect[:, i, s]))  
    
    
    for dof in range(nDOF):
        for s in range(len(xSource)): 
            n0[dof] += (np.dot(dir_vect[:, dof, s], Norm[:, dof]) / 
                    np.linalg.norm(dir_vect[:, dof, s]) / 
                    np.linalg.norm(Norm[:, dof]))
    return n0

def greenMonopole(source, grid, k, position="centroids"):
    """
    Compute Green function between a (NSource, 3) source array and bempp grid object 
    (either centroids or vertices - DOFs)

    Parameters
    ----------
    source : ndarray,
        Position of monopole sources in space.
    grid : bempp-cl grid object,
        Grid containting mesh information.
    k : float
        Wavenumber.
    position : str, optional
        Which DOFs. Can be "centroids" or "vertices". The default is "centroids".

    Returns
    -------
    G : ndarray
        Vector of Green functions.

    """
    if position=="centroids":
        X_grid = grid.centroids
    else:
        X_grid = grid.vertices
    G = np.zeros(len(X_grid), dtype=complex)
    for i in range(len(source)):
        dist = ((source[i, 0]-X_grid[:, 0])**2 + 
                (source[i, 1]-X_grid[:, 1])**2 +
                (source[i, 2]-X_grid[:, 2])**2)**0.5
        G += np.exp(1j*k*dist) / (4*np.pi*dist)
    return G

def buildGridFunction(space, *args):
    dof = space.grid_dof_count
    coeff = np.ones(dof, dtype=complex)
    for arg in args:
        coeff *= arg
    grid_fun = bempp.api.GridFunction(space, coefficients=coeff)
    return grid_fun

def element2vertice_pressure(grid, p_mesh):
    """
    "Transpose" element pressure to vertex pressure.

    Parameters
    ----------
    grid : bempp-cl grid object
        Simulation grid.
    p_mesh : ndarray of gridfunctions
        Result of the pressure over the system's mesh.

    Returns
    -------
    vertice_pressure : ndarray of gridfunctions
        Transposed element to vertex pressure.
    """

    vertices = grid.vertices
    elements = grid.elements
    N = vertices.shape[1]
    M = elements.shape[1]
    
    # Step 1: Initialize an array to accumulate pressure per vertex
    vertice_pressure = np.empty(p_mesh.shape, dtype="object")
    space_vert = bempp.api.function_space(grid, "P", 1)
    
    # Step 2: Loop through each element and distribute pressure
    for f in range(p_mesh.shape[0]):
        for rs in range(p_mesh.shape[1]):
            coeff = np.zeros((len(p_mesh), N), dtype=complex)
            vertex_count = np.zeros(N)  # To count contributions per vertex
            for i in range(M):
                for j in range(3):  # Each element has 3 vertices
                    v_idx = elements[j, i]  # Get vertex index
                    coeff[f, v_idx] += p_mesh[f, rs].coefficients[i]  # sum elements to a single vertex position
                    vertex_count[v_idx] += 1  # Count contributions
        
            nonzero_mask = vertex_count > 0
            coeff[f, nonzero_mask] /= vertex_count[nonzero_mask] 
            vertice_pressure[f, rs] = bempp.api.GridFunction(space_vert, 
                                                             coefficients=coeff[f, :])
    
    # for f in range(p_mesh.shape[0]):
    #     for rs in range(p_mesh.shape[1]):
    #         projection_coeff = p_mesh[f, rs].projections(space_vert)
    #         vertice_pressure[f, rs] = bempp.api.GridFunction(space_vert, coefficients=projection_coeff)
    return vertice_pressure

 
#%% Mesh-related
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

#%% Admittance
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
    absSurf_in = np.array(absorbingSurface)
    surfImp_in = surfaceImpedance
    Nsurf      = len(absSurf_in)             # number of absorbing surfaces
    grid       = spaceP.grid
    dofCount   = spaceP.grid_dof_count

    if absSurf_in.shape[0] == 0 :
        admittanceMatrix = None
    else:
        admittanceMatrix = np.zeros([dofCount, Nfft], dtype=complex) #np.ones([dofCount, Nfft], dtype=complex) * 2.5e-4 # corresponds to 0.05% damping
        for surf in range(Nsurf):
            tmp_surf = absSurf_in[surf]  # current surface on which we apply admittance coefficients
            vertex, _ = get_group_points(grid, tmp_surf)
            for f in range(Nfft):
                try:
                    Yn = rho_0 * c_0 / surfImp_in[surf][f] # rho_0 * c_0
                except:
                    Yn = rho_0 * c_0 / surfImp_in[surf] # rho_0 * c_0
                admittanceMatrix[vertex, f] = np.ones(len(vertex)) * Yn
    return admittanceMatrix 



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



