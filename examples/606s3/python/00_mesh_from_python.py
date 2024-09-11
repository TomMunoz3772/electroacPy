"""
ElectroacPy - 2024-04-24_DEV branch

18/07/2024

This script shows how to load a *.step geometry file and mesh it from Python.
Will need gmsh api, which should already be installed in DEV branch.

The meshing class is set in the generalToolbox module.
"""

# import generalToolbox as gtb      # either load entire toolbox
from generalToolbox import meshCAD  # or just the meshCAD class

# MESH A FULL STEP FILE
path2step = "../geo/step_files/606_PORT.step"
cad = meshCAD(path2step, maxSize=343/2e3/6, minSize=343/2e3/60)     # by default maxSize=0.057 and minSize=maxSize/10
cad.addSurfaceGroup("midrange", surfaces=[8], groupNumber=1)  # surface index is found using a CAD software
cad.addSurfaceGroup("port", surfaces=[7], groupNumber=2)
cad.mesh("../geo/mesh_files/606_port")  # will mesh into corresponding file

# MESH SELECTED FACES ONLY
path2step = "../geo/step_files/606_PORT.step"
cad = meshCAD(path2step)
cad.addSurfaceGroup("midrange", [8], 1)
cad.addSurfaceGroup("baffle", [3], 2)
cad.mesh("../geo/mesh_files/606_baffle", excludeRemaining=True)  # will only mesh "midrange" and "baffle" (surf 8, 3)

# REFINE SURFACE
path2step = "../geo/step_files/606_PORT.step"
cad = meshCAD(path2step)
cad.addSurfaceGroup("midrange", [8], 1, meshSize=343/2500/6)
cad.addSurfaceGroup("baffle", [3], 2)
cad.mesh("../geo/mesh_files/606_baffle_refined", excludeRemaining=True)

"""
Go and check the geo folder for *.msh files, these can be opened in Gmsh and Paraview.
"""