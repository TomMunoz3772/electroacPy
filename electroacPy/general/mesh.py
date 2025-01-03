import gmsh
import numpy as np
from copy import copy

class CAD:
    def __init__(self, file, minSize=0.0057, maxSize=0.057, scaling=0.001, meshAlgo=6):
        """
        Import a geometry to be meshed

        Parameters
        ----------
        file : str
            path to .step file.
        minSize : float, optional
            Minimum mesh size (in m). The default is 0.0057.
        maxSize : float, optional
            Maximum mesh size (in m). The default is 0.057.
        scaling : float, optional
            Scaling of geometry (Gmsh tends to get rid of units, hence mm end up in m). The default is 0.001.
        meshAlgo : int, optional
            Type of mesh. The default is 6 ().

        Returns
        -------
        None.

        """
        
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Geometry.OCCScaling", scaling)
        gmsh.option.setNumber("Mesh.Algorithm", meshAlgo)
        gmsh.option.setNumber("Mesh.MeshSizeMin", minSize)
        gmsh.option.setNumber("Mesh.MeshSizeMax", maxSize)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.QualityType", 2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

        gmsh.merge(file)
        n = gmsh.model.getDimension()
        s = gmsh.model.getEntities(n)
        l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])

        self.file = file
        self.minSize  = minSize
        self.maxSize  = maxSize
        self.scaling  = scaling
        self.meshAlgo = 6
        self.entities = gmsh.model.get_entities()
        self.surface_list = get_tags_by_dimension(self.entities, 2)
        # self.ungrouped_surface = copy(self.surface_list)
        self.physical_groups = []



    def addSurfaceGroup(self, name, surface, groupNumber, meshSize=None):
        if meshSize == None:
            meshSize = self.maxSize
        else:
            meshSize = meshSize
        SURF_GROUP = {"name": name, "surface": surface, "groupNumber": groupNumber, 
                      "meshSize": meshSize, "id": "surface"}
        self.physical_groups.append(SURF_GROUP)
        return None
    
    def addVolumeGroup(self, name, volume, groupNumber, meshSize=None):
        if meshSize == None:
            meshSize = self.maxSize
        else:
            meshSize = meshSize
        VOL_GROUP = {"name": name, "volume": volume, "groupNumber": groupNumber, 
                      "meshSize": meshSize, "id":"volume"}
        self.physical_groups.append(VOL_GROUP)
        return None

    def mesh(self, filename, order=2, excludeRemaining=False, reverseNormals=False):
        # initialization of GMSH

        surface = []
        evalField = np.arange(1, len(self.physical_groups)*2, 2)
        sizeField = evalField + 1
        for i in range(len(self.physical_groups)):
            # group surfaces
            if  self.physical_groups[i]["id"] == "surface":
                surface_tmp = self.physical_groups[i]["surface"]
                groupNumber_tmp = self.physical_groups[i]["groupNumber"]
                name_tmp =  self.physical_groups[i]["name"]
                gmsh.model.geo.add_physical_group(2, surface_tmp, groupNumber_tmp, name_tmp)
                surface.append(surface_tmp)
            elif self.physical_groups[i]["id"] == "volume":
                volume_tmp = self.physical_groups[i]["volume"]
                groupNumber_tmp = self.physical_groups[i]["groupNumber"]
                name_tmp =  self.physical_groups[i]["name"]
                gmsh.model.geo.add_physical_group(3, volume_tmp,
                                                  groupNumber_tmp, name_tmp)
            
        self.ungrouped_surface = filter_surface_list(self.surface_list, surface)
        if excludeRemaining is False:
            gmsh.model.geo.add_physical_group(2, self.ungrouped_surface, -1, "enclosure")
        else:
            pass

        for i in range(len(self.physical_groups)):
            # define mesh field
            # TODO find a way to fix the size modif for volumes
            if self.physical_groups[i]["id"] == "surface":
                gmsh.model.mesh.field.add("MathEval", evalField[i])
                gmsh.model.mesh.field.set_string(evalField[i], "F", 
                                                 str(self.physical_groups[i]["meshSize"]))
                gmsh.model.mesh.field.add("Restrict", sizeField[i])
                gmsh.model.mesh.field.setNumbers(sizeField[i], "SurfacesList", 
                                                 self.physical_groups[i]["surface"])
                gmsh.model.mesh.field.setNumber(sizeField[i], "InField",evalField[i])
            else:
                pass

        gmsh.model.mesh.field.add("Min", evalField[-1]+2)
        gmsh.model.mesh.field.setNumbers(evalField[-1]+2, "FieldsList", list(sizeField))
        gmsh.model.mesh.field.setAsBackgroundMesh(evalField[-1]+2)
        
        
        # ====================================================================
        # extrude PML layer of thickness d_PML and Num_layers layers
        # gmsh.option.setNumber('Geometry.ExtrudeReturnLateralEntities', 0)
        # e = gmsh.model.geo.extrudeBoundaryLayer([(2, 11)], 
        #                                         [5], [0.025],
        #                                         False)

        # # retrieve tags from PML interface (bottom surf), PML subvolumes
        # bottom_surf  = [s[1] for s in [(2, 11)]]
        # pml_entities = [s for s in e if s[0] == 3]
        # pml_volumes  = [s[1] for s in pml_entities]
        
        
        # # set PML volume and PML interface physical groups, generate the entire 3D mesh
        # # and export it in .msh file
        # gmsh.model.addPhysicalGroup(3, pml_volumes, 20, "PML")
        # print("bottom_surf: ", bottom_surf)
        # gmsh.model.addPhysicalGroup(2, bottom_surf, 21, "pml_int")
        # ====================================================================
        
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(order)

        if reverseNormals is True:
            # Get elements and nodes of the mesh
            gmsh.model.mesh.reverse()


        gmsh.write(filename + ".msh")
        gmsh.finalize()
        return None

    def export2stl(self, fileName):
        gmsh.write(fileName + ".stl")
        return None

    def mesh2tmp(self, fileName):
        for i in range(len(self.surface_list)):
            self.addSurfaceGroup(self.surface_list[i], [self.surface_list[i]], self.surface_list[i])
        gmsh.model.mesh.generate(2)
        gmsh.write(fileName + ".msh")
        return None

    def finalize(self):
        gmsh.finalize()


# Useful function
def get_tags_by_dimension(tag_list, target_dimension):
    tags = [tag for dim, tag in tag_list if dim == target_dimension]
    return tags

def filter_surface_list(surface_list, list_of_small_lists):
    exclusion_set = set().union(*list_of_small_lists)
    filtered_list = [item for item in surface_list if item not in exclusion_set]
    return filtered_list