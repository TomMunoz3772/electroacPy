## Boundary Element Modeler
ElectroacPy BEM modeler aims at automating the creation of [bempp-cl studies](https://bempp.com/installation.html): applying velocity to radiating elements, placement of evaluation points in 3D space, definition of the system based on boundary conditions, etc.

The BEM modeler is mostly focused on exterior acoustic radiation (e.g. a loudspeaker radiating in free-air), hence, a lot of functionality that you would find in commercial software are not available, such as interfaces between domains[^interface_note] or structural analysis. However, it is great for:

[^interface_note]: It should be possible to implement it in bempp-cl.

- diffraction studies, 
- signal processing applied to loudspeaker array (beamforming, spherical harmonics, etc.),
- projection of measured/simulated acceleration data.

```{figure} ./drawings/BEM.svg
    :name: 

Structure of the Boundary-Element Modeler.
```

### `.acoustic_study[]` 
The first computational step is to estimate the acoustic pressure on boundaries. The method `.study_acousticBEM()` is used for that effect. The inputs of `.study_acousticBEM()` are:

- `name`, study name, for reference purpose,
- `meshPath`, the path to simulated mesh,
- `acoustic_radiator`, the reference to **enclosure** and / or **driver** object that have been created beforehand, 
- `domain`, relates to the type of study (i.e. `"interior"` or `"exterior"`), default is set to `"exterior"`.

A few `**kwargs` are also available:

- `tol`, tolerance of the GMRES solver,
- `boundary_conditions`, a **boundaryCondition** object which defines infinite boundaries and surfaces impedance,
- `direction`, list of vector that add specific direction coefficients to the radiating surfaces, for example: `[[0, 1, 0]]` for a single driver radiating toward *+y*, or `[[1, 0, 0], False, [1, 0, 0]]` for three drivers, with two radiating toward *+x* and one with normal radiation direction.

In the following code, we define a free-field and a ground-reflected study.

```python
import electroacPy as ep

#%% load data
system = ep.load("03_LEM_Recap")

#%% Define free-field study
system.study_acousticBEM("free-field", 
                         "../geo/mesh/studio_monitor.msh", 
                         ["ported_LF", "TW29", "sealed_MF"], 
                         domain="exterior")

#%% Define boundary conditions
from electroacPy.acousticSim.bem import boundaryConditions

bc = boundaryConditions()
bc.addInfiniteBoundary(normal="z", offset=-1)
system.study_acousticBEM("inf_ground", 
                        "../geo/mesh/studio_monitor.msh",
                        ["ported_LF", "TW29", "sealed_MF"], domain="exterior", 
                        boundary_conditions=bc)

#%% run boundary operators
system.run()

#%% save state
ep.save("04_BEM_setup", system)
```
 
The tweeter pressure on system's boundaries is shown in {numref}`boundary-pressure`. The mirrored mesh is due to the infinite boundary condition on the *xy* plane.

```{figure} ./boundary_images/mesh_pressure_tweeter_b.png
    :name: boundary-pressure

Estimated pressure on boundary.
```

### `.evaluation[]`
While computing the boundary pressure is the most computationally expensive step, it doesn't provide a complete picture of the system. Information about the system's radiated pressure (e.g. directivity, baffle diffraction, etc.) is obtained through the **evaluation** class, which efficiently automates the placement of observation points and visualization of results.

When an acoustic study is created, an **evaluation** object is automatically created. Hence, in regards to evaluations, the "only" functionality of **loudspeakerSystem** is to populate `.evaluation["reference_study"]` with observation setups. For now, there are five different observation types:

- `.evaluation_polarRadiation()`, get directivity information,
- `.evaluation_pressureField()`, pressure across rectangular screens,
- `.evaluation_fieldPoint()`, pressure at one (or more) user defined point,
- `.evaluation_sphericalRadiation()`, acoustic radiation within a sphere of radius $r$,
- `.boundingBox()`, pressure within a parallelepiped.

For our system, we define two polar radiations and two pressure-fields.
```python
import electroacPy as ep

#%% load data
system = ep.load("04_BEM_setup")

#%% Define evaluations
system.evaluation_polarRadiation(["free-field", "inf_ground"], "polar_hor", 
                                 -180, 180, 5, on_axis="x", direction="y",
                                 radius=2, offset=[0, 0, 0.193])

system.evaluation_pressureField(["free-field", "inf_ground"], "field_ver", 
                                L1=3, L2=2, step=343/2500/6, 
                                plane="xz", offset=[-1.5, 0, -1])

system.evaluation_pressureField(["free-field", "inf_ground"], "field_hor", 
                                L1=3, L2=2, step=343/2500/6, 
                                plane="xy", offset=[-1.5, -1, 0.193])

system.evaluation_polarRadiation("free-field", "polar_ver", 
                                -180, 180, 5,  "x", "z", 
                                radius=2, offset=[0, 0, 0.193])

system.evaluation_polarRadiation("inf_ground", "polar_ver", 
                                0, 180, 5, "x", "z",
                                radius=2, offset=[0, 0, 0.193])

# system.plot_system("free-field")
# system.plot_system("inf_ground")

#%% run potential operators and plot results
system.run()
system.plot_results()

#%% save state
ep.save("05_evaluation_setup", system)
```

As you can see in this code, we use the `.plot_system()` function^[Similar to Tkinter window, PyVista plotter will *also* stop execution of code.] --- this will display the 3D placement of evaluation points and boundaries as shown by {numref}`plot-system`. You may also have noticed that the `.run()` command is used again: electroacPy automatically skips any boundary and potential evaluations already computed. 

```{figure} ./boundary_images/system_floor_with_eval_b.png
    :name: plot-system

System under study and related evaluations.
```

Acoustic results are displayed with `.plot_results()`, called from the system object. By default, all results will pop-up. It is possible to return specific results using the following arguments:

- `study`, select specific study or group of study,
- `evaluation`, select specific evaluation(s)
- `radiatingElement`, isolate radiating surfaces / elements
- `bypass_xover`, bypasses filters if **crossovers** are defined.

The following code displays {numref}`plot-field`.
```python
system.plot_results(study="inf_ground", 
                    evaluation=["field_hor", "field_ver"],
                    radiatingElement=[1, 2])
```

```{figure} ./boundary_images/field_plotter_2_b.png
    :name: plot-field

Pressure-Field plotter.
```

When plotting a pressure field, a left click will return the frequency response function of the closest points. {numref}`plot-field-point` is an example of this tool. The legend should display the (x, y, z) position of field point.

```{figure} ./boundary_images/field_point_from_field_b.svg
    :name: plot-field-point

Extraction of FRF from field point.
```

Polar responses are displayed using a Matplotlib viewer. It helps navigating through frequencies / angle with four sub-figures: SPL and normalized directivity, pressure and polar response --- see {numref}`plot-dir`.

```{figure} ./boundary_images/directivity_plotter_b.svg
    :name: plot-dir

Directivity plotter. Vertical radiation in free-field.
```

<!-- In that directivity plot, you can clearly see the limits of the simulation mesh. Because the mesh is relatively coarse, the accuracy of the results decreases significantly at high frequencies, to the point where the off-axis (±180°) pressure is no longer correct. Using a finer mesh will give better results, but will be longer to compute.  -->


### Code Summary
```python
import electroacPy as ep
from electroacPy.acousticSim.bem import boundaryConditions

#%% load data
system = ep.load("03_LEM_Recap")

#%% Define free-field study
system.study_acousticBEM("free-field", 
                         "../geo/mesh/studio_monitor.msh",
                         ["ported_LF", "TW29", "sealed_MF"], domain="exterior")

#%% Define boundary conditions
bc = boundaryConditions()
bc.addInfiniteBoundary(normal="z", offset=-1)
system.study_acousticBEM("inf_ground", 
                        "../geo/mesh/studio_monitor.msh",
                        ["ported_LF", "TW29", "sealed_MF"], domain="exterior", 
                        boundary_conditions=bc)

#%% Define evaluations
system.evaluation_polarRadiation(["free-field", "inf_ground"], "polar_hor", 
                                 -180, 180, 5, "x", "y",
                                 radius=2, offset=[0, 0, 0.193])
system.evaluation_pressureField(["free-field", "inf_ground"], "field_ver", 
                                3, 2, 343/2500/6, "xz", offset=[-1.5, 0, -1])
system.evaluation_pressureField(["free-field", "inf_ground"], "field_hor", 
                                3, 2, 343/2500/6, "xy", 
                                offset=[-1.5, -1, 0.193])
system.evaluation_polarRadiation("free-field", "polar_ver", 
                                -180, 180, 5,  "x", "z", 
                                radius=2, offset=[0, 0, 0.193])
system.evaluation_polarRadiation("inf_ground", "polar_ver", 
                                0, 180, 5, "x", "z",
                                radius=2, offset=[0, 0, 0.193])

#%% run potential operators and plot results
system.run()
system.plot_results()

#%% save state
ep.save("06_acoustic_radiation", system)
```



