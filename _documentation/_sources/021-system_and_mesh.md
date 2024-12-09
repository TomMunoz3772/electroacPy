## System overview and meshing

In the following sections, a three-way monitor is used as an example. {numref}`cad-model` shows the system in its relatively complex CAD form: drivers are fully detailed, including the cone, basket, motor, etc., and the enclosure accounts for internal separations, port dimensions, and wall thicknesses to ensure everything fits properly. However, for our specific case study, this level of detail is unnecessary, as many features will not significantly impact acoustic radiation. Therefore, we simplify the system by retaining only the outer shell and representing each radiating element as flat surfaces --- see {numref}`sim-model`. These significant simplifications should still provide valuable insights into acoustic radiation.


```{figure} ./system_images/studio_monitor.png
    :name: cad-model
    :width: 300px
    :align: left

CAD model.
```

```{figure} ./system_images/studio_monitor_sim.png
    :name: sim-model
    :width: 300px
    :align: right

Simplified model for simulations.
```

### Mesh
A mesh of the system is required for simulation with bempp. As noted in the bempp documentation, mesh imports are handled with [meshio](https://github.com/nschloe/meshio). Compatible mesh formats are listed on its GitHub page. To simplify the process of meshing geometries, ElectroacPy provides a wrapper for the [Gmsh](https://gmsh.info/) API. This is what we'll use for our study.

Once your geometry is exported as a `.step` file, you can create a `CAD` object with `generalToolbox`:

```python
import generalToolbox as gtb

cad = gtb.meshCAD("../geo/step_export/simulation_cad.step")
```
By default, the maximum mesh size corresponds to approximately 6 elements per wavelength at a frequency of 1 kHz, which equals 57 mm for a sound speed of 343 m/s. The minimum mesh size is ten times smaller: 5.7 mm. You can change these values by passing the `minSize` and `maxSize` arguments when importing a geometry:

```python
cad = gtb.meshCAD("../geo/step_export/simulation_cad.step", minSize=343/2500/60, 
maxSize=343/2500/6)
```

Next, surfaces must be grouped and referenced to set independent boundary or radiation conditions on separate surfaces. In our example, four radiating surfaces are considered:


```python
cad.addSurfaceGroup("woofer", surface=[8], groupNumber=1)
cad.addSurfaceGroup("port", [10], 2)
cad.addSurfaceGroup("midrange", [9], 3)
cad.addSurfaceGroup("tweeter", [7], 4)
```

This code creates four surface groups numbered from `1` to `4`, each referencing a different surface in the geometry. Surface indices can be retrieved from your CAD software. For example, FreeCAD displays the element index in the bottom-left corner when hovering over a surface, line, or point.

```{figure} ./system_images/pre_mesh.png
    :name: select-surface

Getting surface index.
```

Once grouping of surfaces is done, run the `mesh` command. By default, it creates a 2D mesh under Gmsh's `.msh` file format. 
```python
cad.mesh("../geo/mesh/studio_monitor")
```

```{figure} ./system_images/mesh_1kHz.png
    :name: mesh-1kHz

Geometry meshed with a 1 kHz resolution.
```

Of course, using a global mesh size can be restrictive. Therefore, it is possible to specify separate mesh sizes. In the following code, the global mesh size is set for 1 kHz, while the mesh is refined to 5 kHz over the tweeter and baffle surfaces.

```python
#%% Set global size (1kHz)
lmax = 343/1e3/6
lmin = lmax/10
c
#%% import geo
cad = gtb.meshCAD("../geo/step_export/simulation_cad.step", 
                  minSize=lmin, maxSize=lmax)
                  
#%% group surfaces
cad.addSurfaceGroup("woofer", [8], 1)
cad.addSurfaceGroup("port", [10], 2)
cad.addSurfaceGroup("midrange", [9], 3)
cad.addSurfaceGroup("tweeter", [7], 4, meshSize=343/5e3/6)
cad.addSurfaceGroup("baffle", [3], 5, meshSize=343/5e3/6)
cad.mesh("../geo/mesh/studio_monitor_refined")
```

```{figure} ./system_images/mesh_refined.png
    :name: mesh-refined

Refined mesh.
```

As you can see, although the smaller mesh was intended to be applied only to the tweeter and baffle, it is also present on the subwoofer, midrange, and port. This is a potential area for improvement in ElectroacPy's automated API calls.

A final note: after creating surface groups, `meshCAD` will automatically regroup all remaining surfaces into a single group called "enclosure" --- this ensures that all surfaces are properly loaded into bempp.
